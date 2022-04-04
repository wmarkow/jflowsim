//******************************************************************
//                                                                 
//  DistillerBuilder.java                                               
//  Copyright 2022 PSI AG. All rights reserved.              
//  PSI PROPRIETARY/CONFIDENTIAL. Use is subject to license terms
//                                                                 
// ******************************************************************

package jflowsim.controller.solverbuilder;

import jflowsim.model.numerics.Solver;
import jflowsim.model.numerics.UniformGrid;
import jflowsim.model.numerics.lbm.distiller.DistillerSolver;
import jflowsim.model.numerics.lbm.testcases.TestCase;
import jflowsim.model.numerics.lbm.testcases.distiller.BreakingDamTestCase;
import jflowsim.model.numerics.lbm.testcases.distiller.DistillerSmallDropTestCase;
import jflowsim.model.numerics.lbm.testcases.distiller.DistillerTestCase;

public class DistillerBuilder extends SolverBuilder {

    public DistillerBuilder() {
	testCaseSet.put(DistillerTestCase.class.getSimpleName(), new DistillerTestCase());
	testCaseSet.put(DistillerSmallDropTestCase.class.getSimpleName(), new DistillerSmallDropTestCase());
	testCaseSet.put(BreakingDamTestCase.class.getSimpleName(), new BreakingDamTestCase());
    }

    @Override
    public Solver createSolver(UniformGrid grid) {
	return new DistillerSolver(grid);
    }

    @Override
    public UniformGrid createGrid(String testcase) {
	TestCase testCaseCreator = testCaseSet.get(testcase);

	return testCaseCreator.getGrid();
    }

}
