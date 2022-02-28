package jflowsim.model.numerics;

import java.util.ArrayList;

import jflowsim.model.ModelManager;
import jflowsim.model.geometry2d.Geometry2D;
import jflowsim.model.geometry2d.utilities.Rounding;
import jflowsim.model.numerics.utilities.GridNodeType;
import jflowsim.view.headupdisplay.HeadUpDisplay;

/***
 * An abstract implementation of lattice grid.
 */
public abstract class UniformGrid {

    /***
     * Number of x and y discrete cells along X and Y axis.
     */
    public int nx, ny;

    /***
     * Physical bounds of the grid.
     */
    private double minX, minY, maxX, maxY;

    /***
     * Enables/disables periodic boundary conditions for X and Y axes. Example: if
     * periodicX is true then the fluid at the east side of the grid may disappear,
     * and appear on the west side of the grid.
     */
    public boolean periodicX = false, periodicY = false;

    /***
     * Physical dimensions of the grid.
     */
    private double length, width; // PRIVATE, damit Zugriff nur ueber getter/setter zwecks min/max-update

    /***
     * dx is a physical size (width and height) of the grid cell. Grid cell is a
     * square. dt is a physical smallest amount of time used in the simulation dv is
     * the smallest velocity ?
     */
    public double dx, dt, dv;

    /***
     * Viscosity of the grid. TODO: why is this a property of the whole grid and not
     * inside of each cell?
     */
    public double viscosity;

    /***
     * Gravity along X and Y axes.
     */
    public double gravityX, gravityY;

    /***
     * A general gravity value inside of the grid.
     */
    public double gravity = 9.81;

    /***
     * Used as a hint on how often update the statistics information (i.e. by
     * calling {@link Solver#update()}. This is expressed in a discrete time step
     * unit.
     */
    public int updateInterval = 25;

    /***
     * Discrete time step used during simulation.
     */
    public int timestep;

    /***
     * Real time in seconds used during simulation.
     */
    public double real_time; // real in secs

    /***
     * MNU per seconds? What is this?
     */
    public double mnups;

    /***
     * A name of the simulation test case.
     */
    public String testcase;

    /***
     * A list of boundary conditions.
     */
    public ArrayList<BoundaryCondition> bcList = new ArrayList<BoundaryCondition>();

    /***
     * A one dimensional table which stores the specific type of cell. The index of
     * cell at coordinates (x,y) = y * nx + x
     */
    public int type[];

    public UniformGrid(double _length, double _width, double _dx) {

	this.setLength(_length);
	this.setWidth(_width);
	this.dx = _dx;

	this.nx = (int) Math.round(this.getLength() / this.dx) + 1;
	this.ny = (int) Math.round(this.getWidth() / this.dx) + 1;

	this.allocateMemory();
    }

    public UniformGrid(double _length, double _width, int _nx, int _ny) {

	this.setLength(_length);
	this.setWidth(_width);

	this.nx = _nx;
	this.ny = _ny;

	this.dx = this.getLength() / (this.nx - 1);

	this.allocateMemory();
    }

    protected abstract void allocateMemory();

    public abstract void updateParameters();

    public abstract void refineGrid(double scaleFactor);

    public void map2Grid(ModelManager model) {

	for (int x = 0; x < nx; x++) {
	    for (int y = 0; y < ny; y++) {
		if (getType(x, y) == GridNodeType.SOLID) {
		    setType(x, y, GridNodeType.FLUID);
		}
	    }
	}
	// map Geometry to Grid
	for (Geometry2D geo : model.getGeometryList()) {
	    geo.map2Grid(this);
	}
    }

    public void addBC(BoundaryCondition bc) {
	bcList.add(bc);
    }

    public void updateHeadUpDisplay(HeadUpDisplay hud) {
    }

    public abstract double getScalar(int x, int y, int scalar_name);

    // checks if (x,y) is inside the grid
    public boolean isPointInside(double x, double y) {
	if (x < minX || x > maxX) {
	    return false;
	}
	if (y < minY || y > maxY) {
	    return false;
	}
	return true;
    }

    public void setLength(double length) {
	this.minX = 0.0;
	this.maxX = length;

	this.length = length;
    }

    public void setGridResolution(double deltaX) {
	this.dx = deltaX;
    }

    public void setWidth(double width) {
	this.minY = 0.0;
	this.maxY = width;

	this.width = width;
    }

    public double getLength() {
	return this.length;
    }

    public double getWidth() {
	return this.width;
    }

    public double getMinX() {
	return minX;
    }

    public double getMaxX() {
	return maxX;
    }

    public double getMinY() {
	return minY;
    }

    public double getMaxY() {
	return maxY;
    }

    public int transCoord2XIndex(double x, int mode) {
	int ix = Rounding.round(x / dx + 1.0, mode);
	if (ix < 0) {
	    ix = 0;
	}
	if (ix > nx - 1) {
	    ix = nx - 1;
	}
	return ix;
    }

    public int transCoord2YIndex(double y, int mode) {
	int iy = Rounding.round(y / dx + 1.0, mode);
	if (iy < 0) {
	    iy = 0;
	}
	if (iy > ny - 1) {
	    iy = ny - 1;
	}
	return iy;
    }

    public double transXIndex2Coord(int x) {
	return minX + dx * x;
    }

    public double transYIndex2Coord(int y) {
	return minY + dx * y;
    }

    public int getType(int x, int y) {
	return this.type[y * this.nx + x];
    }

    public void setType(int x, int y, int type) {
	this.type[y * this.nx + x] = type;
    }

    public void setViscosity(double _viscosity) {
	this.viscosity = _viscosity;
    }

    public void setGravity(double _gravityX, double _gravityY) {
	this.gravityX = _gravityX;
	this.gravityY = _gravityY;
    }

    public void setTimeStep(double _dt) {
	this.dt = _dt;
    }
}
