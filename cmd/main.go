package main

import (
	"fmt"
	"image"
	"image/color"
	"log"
	"math"
	"runtime"
	"sync"

	"github.com/hajimehoshi/ebiten/v2"
	"github.com/hajimehoshi/ebiten/v2/ebitenutil"
	"github.com/hajimehoshi/ebiten/v2/inpututil"
)

const (
	screenWidth   = 800
	screenHeight  = 600
	maxIterations = 256 // Default maximum iterations for coloring
)

// Game implements the ebiten.Game interface.
type Game struct {
	// Fractal State
	Center  complex128 // Current center point (C_real + C_imag*i)
	Zoom    float64    // Current zoom level (pixels per complex unit)
	MaxIter int        // Current maximum iterations

	// Rendering State
	MandelbrotImage *ebiten.Image // The cached image of the fractal
	rawImage        *image.RGBA   // The raw RGBA image data being calculated

	// Concurrency and Synchronization
	mu                 sync.Mutex // Mutex to protect shared data (like rawImage)
	calculationRunning bool       // Flag to prevent multiple concurrent calculations
	needsRecalculation bool       // Flag set by Update when input changes the view

	// Input State
	dragStartPos image.Point // Starting position of a mouse drag
	isDragging   bool        // Whether the mouse is currently being dragged
}

func NewGame() *Game {
	g := &Game{
		Center:             complex(-0.5, 0), // Start at the main cardiod center
		Zoom:               200.0,            // Initial zoom (200 pixels cover 1 unit in complex space)
		MaxIter:            maxIterations,
		rawImage:           image.NewRGBA(image.Rect(0, 0, screenWidth, screenHeight)),
		MandelbrotImage:    ebiten.NewImage(screenWidth, screenHeight),
		needsRecalculation: true, // Calculate the first frame immediately
	}
	// Copy the initial empty raw image to the Ebiten image
	g.MandelbrotImage.WritePixels(g.rawImage.Pix)
	return g
}

// handleZoomChange updates the zoom level and max iterations based on a factor,
// keeping the specified complex point 'c' constant. If 'c' is zero, it centers zoom on the screen center.
// This logic is factored out to be used by both mouse scroll and keyboard input.
func (g *Game) handleZoomChange(c complex128, factor float64) {
	// If the complex point is 0, zoom is centered on the current screen center (g.Center)
	if c == 0 {
		g.Zoom *= factor
	} else {
		// Re-center the view around the cursor point 'c' after zooming
		g.Zoom *= factor
		// This keeps the point under the cursor stable during zoom.
		g.Center = g.Center + (g.Center-c)*(complex(factor-1, 0))/complex(factor, 0)
	}

	// Adjust max iterations based on zoom (deep zoom requires more iterations)
	g.MaxIter = int(math.Max(maxIterations, 5*math.Log2(g.Zoom)*20))
	if g.MaxIter > 2000 {
		g.MaxIter = 2000
	} else if g.MaxIter < maxIterations {
		g.MaxIter = maxIterations
	}
}

// Update handles game logic and input.
func (g *Game) Update() error {
	var viewChanged bool

	// --- 1. Handle Zoom (Scroll Wheel) ---
	_, scrollY := ebiten.Wheel()
	if scrollY != 0 {
		cursorX, cursorY := ebiten.CursorPosition()
		// Get the complex point under the cursor
		c := g.screenToComplex(cursorX, cursorY)

		zoomFactor := 1.1 // Zoom by 10% per scroll tick
		if scrollY < 0 {
			zoomFactor = 1.0 / zoomFactor // Zoom out
		}

		g.handleZoomChange(c, zoomFactor)

		viewChanged = true
	}

	// --- 2. Handle Panning (Click and Drag) ---
	if inpututil.IsMouseButtonJustPressed(ebiten.MouseButtonLeft) {
		g.dragStartPos.X, g.dragStartPos.Y = ebiten.CursorPosition()
		g.isDragging = true
	} else if inpututil.IsMouseButtonJustReleased(ebiten.MouseButtonLeft) {
		g.isDragging = false
		viewChanged = true // A drag always ends in a recalculation
	} else if g.isDragging {
		// Dragging in progress: update the center based on mouse movement
		currentX, currentY := ebiten.CursorPosition()

		// Calculate the difference in complex space
		startC := g.screenToComplex(g.dragStartPos.X, g.dragStartPos.Y)
		currentC := g.screenToComplex(currentX, currentY)

		// The new center is the old center minus the displacement vector
		g.Center -= (currentC - startC)

		g.dragStartPos.X = currentX
		g.dragStartPos.Y = currentY
		viewChanged = true
	}

	// --- 3. Handle Keyboard Input for Zoom and Recenter ---
	// 'Z' to Recenter (Resets all fractal state)
	if inpututil.IsKeyJustPressed(ebiten.KeyZ) {
		g.Center = complex(-0.5, 0) // Reset to the default center
		g.Zoom = 200.0              // Reset to the initial zoom
		g.MaxIter = maxIterations
		viewChanged = true
	}

	// '+' (Key Equal) to Zoom In (centered)
	if inpututil.IsKeyJustPressed(ebiten.KeyEqual) {
		// Use complex(0, 0) to indicate zoom should be centered on g.Center
		g.handleZoomChange(0, 1.1)
		viewChanged = true
	}

	// '-' (Key Minus) to Zoom Out (centered)
	if inpututil.IsKeyJustPressed(ebiten.KeyMinus) {
		// Use complex(0, 0) to indicate zoom should be centered on g.Center
		g.handleZoomChange(0, 1.0/1.1)
		viewChanged = true
	}

	// --- 4. Trigger Recalculation ---
	if viewChanged && !g.calculationRunning {
		g.needsRecalculation = true
	}

	if g.needsRecalculation && !g.calculationRunning {
		g.calculationRunning = true
		g.needsRecalculation = false // Calculation is starting

		// Launch the heavy calculation in a new goroutine
		go g.calculateMandelbrotAsync()
	}

	return nil
}

// Draw renders the current game state.
func (g *Game) Draw(screen *ebiten.Image) {
	// Draw the cached Mandelbrot image
	screen.DrawImage(g.MandelbrotImage, &ebiten.DrawImageOptions{})

	// Draw status text for debugging and feedback
	status := fmt.Sprintf(
		"TPS: %0.2f\nFPS: %0.2f\nMax Iter: %d\nCenter: (%0.8f, %0.8f)\nZoom: %0.2f\nStatus: %s",
		ebiten.ActualTPS(),
		ebiten.ActualFPS(),
		g.MaxIter,
		real(g.Center),
		imag(g.Center),
		g.Zoom,
		g.getCalculationStatus(),
	)
	ebitenutil.DebugPrint(screen, status)
}

// Layout is called when the window size changes.
func (g *Game) Layout(outsideWidth, outsideHeight int) (int, int) {
	return screenWidth, screenHeight
}

// screenToComplex converts a pixel coordinate to a complex number.
func (g *Game) screenToComplex(px, py int) complex128 {
	// Calculate complex coordinate based on screen center and zoom
	x := float64(px-screenWidth/2)/g.Zoom + real(g.Center)
	y := float64(py-screenHeight/2)/g.Zoom + imag(g.Center)
	return complex(x, y)
}

// getCalculationStatus returns a human-readable status message.
func (g *Game) getCalculationStatus() string {
	if g.calculationRunning {
		return "Calculating..."
	}
	if g.needsRecalculation {
		return "Ready to Calculate"
	}
	return "Idle"
}

// calculateMandelbrotAsync runs the fractal generation in a background goroutine.
func (g *Game) calculateMandelbrotAsync() {
	defer func() {
		g.mu.Lock()
		g.calculationRunning = false
		g.mu.Unlock()
	}()

	// Create a new image for this calculation round
	newImage := image.NewRGBA(image.Rect(0, 0, screenWidth, screenHeight))

	// Get the current fractal parameters safely
	g.mu.Lock()
	center := g.Center
	zoom := g.Zoom
	maxIter := g.MaxIter
	g.mu.Unlock()

	// Use a WaitGroup for parallel processing
	var wg sync.WaitGroup
	// Optimal Go routine count is often around the number of CPU cores
	numWorkers := runtime.NumCPU()

	// Divide the image into vertical strips for parallel calculation
	stripWidth := screenWidth / numWorkers

	for i := 0; i < numWorkers; i++ {
		wg.Add(1)

		// Define the start and end pixel columns for this worker
		startX := i * stripWidth
		endX := (i + 1) * stripWidth
		if i == numWorkers-1 {
			endX = screenWidth // Ensure the last worker handles any remainder
		}

		go func(startX, endX int) {
			defer wg.Done()

			// Worker calculates its assigned strip
			for py := 0; py < screenHeight; py++ {
				for px := startX; px < endX; px++ {
					// Map pixel coordinate to complex number c
					x := float64(px-screenWidth/2)/zoom + real(center)
					y := float64(py-screenHeight/2)/zoom + imag(center)
					c := complex(x, y)

					// Calculate divergence
					iter, magSq := g.mandelbrot(c, maxIter)

					// Set color in the temporary image
					newImage.Set(px, py, g.colorize(iter, magSq))
				}
			}
		}(startX, endX)
	}

	wg.Wait() // Wait until all workers complete their strips

	// Once calculation is done, update the main image cache
	g.mu.Lock()
	g.rawImage = newImage
	g.MandelbrotImage.WritePixels(newImage.Pix)
	g.mu.Unlock()
}

// mandelbrot implements the core z = z^2 + c iteration.
func (g *Game) mandelbrot(c complex128, maxIter int) (int, float64) {
	var z complex128
	for n := 0; n < maxIter; n++ {
		// z = z^2 + c
		z = z*z + c

		// Check magnitude squared against escape radius (4)
		magSq := real(z)*real(z) + imag(z)*imag(z)
		if magSq > 4 {
			return n, magSq // Diverged at iteration n
		}
	}
	return maxIter, 0.0 // Inside the set
}

// colorize generates a color based on the iteration count and final magnitude.
func (g *Game) colorize(iter int, magSq float64) color.Color {
	if iter == g.MaxIter {
		return color.Black // Inside the set
	}

	// Smooth coloring calculation (to avoid banding)
	// log2(log(r)/log(2)) where r is the final magnitude (sqrt(magSq))
	// log(r) = 0.5 * log(magSq)
	mu := float64(iter) + 1 - math.Log(0.5*math.Log(magSq))/math.Log(2)

	// Map the continuous iteration count (mu) to a color palette (HSL-like)
	// We use a cosine interpolation to generate vibrant, smooth colors.

	// Simple color cycle based on mu
	red := uint8(math.Sin(0.016*mu+0.0)*230 + 25)
	green := uint8(math.Sin(0.016*mu+0.5)*230 + 25)
	blue := uint8(math.Sin(0.016*mu+1.0)*230 + 25)

	return color.RGBA{R: red, G: green, B: blue, A: 255}
}

func main() {

	ebiten.SetWindowSize(screenWidth, screenHeight)
	ebiten.SetWindowTitle("Interactive Mandelbrot Set (Ebitengine)")

	if err := ebiten.RunGame(NewGame()); err != nil {
		log.Fatal(err)
	}
}
