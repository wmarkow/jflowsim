//******************************************************************
//                                                                 
//  DistillerTestCase.java                                               
//  Copyright 2022 PSI AG. All rights reserved.              
//  PSI PROPRIETARY/CONFIDENTIAL. Use is subject to license terms
//                                                                 
// ******************************************************************

package jflowsim.model.numerics.lbm.testcases.distiller;

import jflowsim.model.numerics.UniformGrid;
import jflowsim.model.numerics.lbm.LbEQ;
import jflowsim.model.numerics.lbm.distiller.DistillerGrid;
import jflowsim.model.numerics.lbm.testcases.TestCase;
import jflowsim.model.numerics.utilities.GridNodeType;

public class BreakingDamTestCase extends TestCase {

    @Override
    public UniformGrid getGrid() {
	double dx = 0.002;
	DistillerGrid grid = new DistillerGrid(0.5715 /* length */, 0.142875 /* height */, dx /* dx */);

	grid.testcase = this.getClass().getSimpleName();

	grid.setLBParameters(0.005 /* nue_lbm */, 0.0 /* forcingX */, -0.001 /* forcingY */);

	// Fill the whole grid with GAS
	for (int i = 0; i < grid.nx * grid.ny; i++) {
	    grid.type[i] = GridNodeType.GAS;
	}

	// fill with water
	int startX = 1;
	int startY = 1;
	int endX = (int) (0.05715 / dx);
	int endY = (int) (0.1143 / dx);
	for (int x = startX; x <= endX; x++) {
	    for (int y = startY; y <= endY; y++) {
		grid.setType(x, y, GridNodeType.FLUID);
		grid.setFill(x, y, 1.0);
	    }
	}

	// Create an interface between fluid and gas.
	for (int x = 0; x < grid.nx; x++) {
	    for (int y = 0; y < grid.ny; y++) {
		// Iterate over the whole grid.

		if (grid.getType(x, y) == GridNodeType.FLUID) {
		    // We have a FLUID node...
		    Boolean changeToInterface = false;

		    for (int dir = 0; dir < 9; dir++) {
			if (grid.getType(x + LbEQ.ex[dir], y + LbEQ.ey[dir]) == GridNodeType.GAS) {
			    // ...which has a direct contact to GAS...
			    changeToInterface = true;
			}
		    }

		    if (changeToInterface) {
			// .. change it to interface but keep in mind that a fluid level is tracked
			// under fill level variable
			grid.setType(x, y, GridNodeType.INTERFACE);
			grid.setFill(x, y, 1.0);
		    }
		}
	    }
	}

	// Boundary conditions
	for (int x = 0; x < grid.nx; x++) {
	    // bottom edge of the column
	    grid.type[x] = GridNodeType.BOUNDARY;
	    // top edge of the column
	    grid.type[x + (grid.ny - 1) * grid.nx] = GridNodeType.BOUNDARY;
	}
	for (int y = 0; y < grid.ny; y++) {
	    // left side of the column
	    grid.type[0 + y * grid.nx] = GridNodeType.BOUNDARY;
	    // right side of the column
	    grid.type[grid.nx - 1 + y * grid.nx] = GridNodeType.BOUNDARY;
	}

	grid.viscosity = 0.05;
	grid.dv = grid.viscosity / grid.nue_lbm * grid.nx / grid.getLength();
	grid.dt = grid.dx / grid.dv;

	// Initial conditions. Fill every node with initial BGK state.
	for (int i = 0; i < grid.nx; i++) {
	    for (int j = 0; j < grid.ny; j++) {

		double[] feq = new double[9];

		LbEQ.getBGKEquilibrium(1.0, 0., 0., feq);

		for (int dir = 0; dir < 9; dir++) {
		    grid.f[(i + j * grid.nx) * 9 + dir] = feq[dir];
		    grid.ftemp[(i + j * grid.nx) * 9 + dir] = feq[dir];
		}

	    }
	}

	return grid;
    }

}
