/**
 * 
 */
package muscle;

import java.util.ArrayList;
import java.util.List;

import com.jidesoft.comparator.BooleanComparator;

import repast.simphony.context.Context;
import repast.simphony.engine.schedule.ScheduledMethod;
import repast.simphony.query.space.grid.GridCell;
import repast.simphony.query.space.grid.GridCellNgh;
import repast.simphony.query.space.grid.MooreQuery;
import repast.simphony.random.RandomHelper;
import repast.simphony.space.SpatialMath;
import repast.simphony.space.continuous.ContinuousSpace;
import repast.simphony.space.continuous.NdPoint;
import repast.simphony.space.grid.Grid;
import repast.simphony.space.grid.GridPoint;
import repast.simphony.util.ContextUtils;
import repast.simphony.util.SimUtilities;

/**
 * @author Andy N Vo
 *
 */

public class Macrophage {
	
	private static ContinuousSpace<Object> space;
	private static Grid<Object> grid;
	private int phenotype; // 0 = M0, 1 = M1, 2 = M2, 3 = resident
	private int age;
	private int phagocytosis;
	private static double[] growthFactors;
	private static double alpha;
	private Object buddy;
	private boolean linked;
	private int stay;
	private int m1num;	// number of M1s
	private int m2num;	// number of M2s
	private int resMnum; // number of resMs
	
	public Macrophage(ContinuousSpace<Object> space, Grid<Object> grid, int phenotype, int age, int phagocytosis, boolean linked, Object buddy) {
		this.space = space;
		this.grid = grid;
		this.phenotype = phenotype;
		this.age = age;
		this.phagocytosis = phagocytosis;
		this.linked = linked;
		this.buddy = buddy;
	}


	// MACROPHAGE RECRUITMENT	
	@ScheduledMethod(start = 2, interval = 1, pick = 1)
	public void macRecruitment() {
		
		Context context = ContextUtils.getContext(this);
		growthFactors = GrowthFactors.getGrowthFactors();
		/*
		 0 = tgf, 1 = tnf, 2 = igf1, 3 = pdgf, 4 = mmpx, 5 = colX, 6 = ECM prot, 7 = il1, 8 = il8, 9 = cxcl2, 10 = cxcl1
		 11 = ccl3, 12 = ccl4, 13 = il6, 14 = mcp, 15 = ifn, 16 = lactoferins, 17 = hgf, 18 = vegf, 19 = mmp12, 20 = gcsf
		 21 = il10, 22 = lipoxins, 23 = resolvins, 24 = ccl17, 25 = ccl22, 26 = colIV, 27 = pge2, 28 = ROS, 29 = fgf, 30 = il4
		*/
		double azurocidin = 0;
		double nitricoxide = 0;
		double ll37 = 0;
		double ccl2 = 0;
		double ccl6 = 0;
		double cx3cl1 = 0;
		double Il13 = 0;
		double recruit = azurocidin + ll37 + ccl2 + growthFactors[12] + growthFactors[24] + growthFactors[25]
				+ growthFactors[11] + ccl6 + growthFactors[13] + cx3cl1 + +growthFactors[14];		//Azurocidin + LL37 + CCL2 + CCL4 + CCL17 + CCL22 + CCL3 + CCL6 + sum [IL6] of patches + CX3CL1  + sum [MCP] of patches
		double deter = growthFactors[27] + growthFactors[22] + nitricoxide + growthFactors[0];		//PGE2 + Lipoxins + NO + sum [TGF] of patches
		double macProb = 0;
		double differential = recruit - deter;
		//System.out.println(differential);
		if (differential > 0) { // threshold to recruit macrophage
			macProb = RandomHelper.nextIntFromTo(0, (int) (15 - 8 + 2 * 6)/10) ; // p(recruit Mac) based on the amounts of DAMPs in the area
			//macProb = RandomHelper.nextIntFromTo(0, (int) Math.ceil(differential / 100)) ;
			//macProb = RandomHelper.nextIntFromTo(0, (int) Math.ceil(differential)) ;
		}
		//if (macProb <= 0 && ((recruit * 1.5) > deter)) { // if there isn't a HUGE deter signal
			//macProb = RandomHelper.nextIntFromTo(0, (5 + 6 / 2 - 3)); // some Mac will be recruited
			//macProb = RandomHelper.nextIntFromTo(0, 2); // some Mac will be recruited
		//}
		//System.out.println(macProb);

		double m1Chance = growthFactors[15] + growthFactors[1] / 2 - growthFactors[21] - growthFactors[0];		//IFN + sum [TNF] of patches / 2 - sum [IL10] of patches - sum [TGF] of patches
		if (m1Chance < 0) {
			m1Chance = 0.01;
		}
		double m2Chance = (growthFactors[21] * 2 + growthFactors[30] + Il13 - growthFactors[15]);		//(sum [IL10] of patches * 2) + IL4 + IL13 - IFN
		if (m2Chance < 0) {
			m2Chance = 0.01;
		}
		/*System.out.println("m1Chance: ");
		System.out.println(m1Chance);
		System.out.println("m2Chance: ");
		System.out.println(m2Chance);*/
				
		m1num = getM1(context).size(); // get number of M1s
		m2num = getM2(context).size(); // get number of M2s
		if (m1num > 1000 || m2num > 1000) {
			macProb = 0; // don't recruit Macs if M1 OR M2 > 100 (NetLogo: 1000 instead of 100)
		}
		while (macProb > 0) {
			double chance = RandomHelper.nextDoubleFromTo(0, m1Chance + m2Chance);
			List<Object> ecms = ECM.getECM(context);
			List<Object> necrosis = Necrosis.getNecrosis(context);
			List<Object> patch = new ArrayList<Object>();
			for (Object obj : ecms) {
				if (obj instanceof ECM && ((ECM) obj).getCollagen() > 0) {
					patch.add(obj);
				}
			} 
			for (Object obj : necrosis) {
				if (obj instanceof Necrosis) {
					patch.add(obj);
				}
			}
			
			if (chance < m1Chance) {
				int index = RandomHelper.nextIntFromTo(0, patch.size() - 1);
				Object spawnpt = patch.get(index);
				GridPoint gridpt = grid.getLocation(spawnpt);
				// NdPoint spacept = space.getLocation(spawnpt);
				Macrophage m1mac = new Macrophage(space, grid, 1, 0, 0, false, null); // phenotype = 1; age = 0; phagocytosis = 0; linked = false; buddy = null
				context.add(m1mac);
				grid.moveTo(m1mac, gridpt.getX(), gridpt.getY());
				// space.moveTo(m1mac, spacept.getX(), spacept.getY());
			}
			else {
				int index = RandomHelper.nextIntFromTo(0, patch.size() - 1);
				Object spawnpt = patch.get(index);
				GridPoint gridpt = grid.getLocation(spawnpt);
				// NdPoint spacept = space.getLocation(spawnpt);
				Macrophage m2mac = new Macrophage(space, grid, 2, 0, 0, false, null); // phenotype = 2; age = 0; phagocytosis = 0; linked = false; buddy = null
				context.add(m2mac);
				grid.moveTo(m2mac, gridpt.getX(), gridpt.getY());
				// space.moveTo(m2mac, spacept.getX(), spacept.getY());
			}
			macProb = macProb - 1;
		}
		alpha = deter;
	}

	@ScheduledMethod(start = 2, interval = 1)
	public void step() {
		// Sense environment
		Context context = ContextUtils.getContext(this);
		growthFactors = GrowthFactors.getGrowthFactors();
		/*
		 0 = tgf, 1 = tnf, 2 = igf1, 3 = pdgf, 4 = mmpx, 5 = colX, 6 = ECM prot, 7 = il1, 8 = il8, 9 = cxcl2, 10 = cxcl1
		 11 = ccl3, 12 = ccl4, 13 = il6, 14 = mcp, 15 = ifn, 16 = lactoferins, 17 = hgf, 18 = vegf, 19 = mmp12, 20 = gcsf
		 21 = il10, 22 = lipoxins, 23 = resolvins, 24 = ccl17, 25 = ccl22, 26 = colIV, 27 = pge2, 28 = ROS, 29 = fgf, 30 = il4
		*/
				
		// Age macrophage
		this.age = this.age + 1;
		//System.out.println(this.phagocytosis);
		// Determine behaviors
		if (this.phenotype == 1) {
			if (this.phagocytosis == 0) {		// hasn't eating anything
				m1Behave(context);
				//System.out.println(this.phagocytosis);
			}
			else if (this.phagocytosis > 0) {		// has eaten apoptotic neutrophil
				m1ApopEat(context);
			}
			else if (this.phagocytosis < 0) {
				//m1DebrisEat(context);
			}
		}
		else if (this.phenotype == 2) {
			m2Behave(context);
		}
		else if (this.phenotype == 3) {
			//resMacBehave(context);
			int dRMdt = (int)((.1 * InflamCell.rmBasal * .3 - .1 * getMres(context).size() * 1.8)*-10);
			int random = RandomHelper.nextIntFromTo(0, dRMdt);
			//System.out.print(random);
			if (random==1) {
					macDie();
				}
		}
		
		if (this.phenotype != 3) {
			// Determine proliferation
			int random = RandomHelper.nextIntFromTo(0, 43);	// originally, NetLogo: random 10 + 5 * 6.7
			m1num = getM1(context).size();	// get number of M1s
			m2num = getM2(context).size();	// get number of M2s
			resMnum = getMres(context).size(); // get number of resMs
			int macCount = m1num + m2num + resMnum; // total macrophages
			if ((random == 1)) {//&& (macCount < 100)) { // NetLogo: 1000 instead of 100
				//macProliferate();
			}
			
			// Determine death
			if (this.age > 10 + 5 * 6.2) {
				if (RandomHelper.nextIntFromTo(0, 4) == 1) { // originally, NetLogo: random 5
					//macDie(context);
					macDie();
				}
			}
		}
	}
	
	// Resident Mac Behavior
	public void resMacBehave(Context<Object> context) {
		//Context context = ContextUtils.getContext(this);
		// if you see DAMPS or IL1 - start producing Cytokines - omit here
		// find spots with collagen
		List<Object> ecmElems = ECM.getECM(context);
		List<Object> ecmWithCol = new ArrayList<Object>();
		for (Object obj : ecmElems) {
			if(((ECM) obj).getCollagen() > 0) {
				ecmWithCol.add(obj);
			}
		}
		if (RandomHelper.nextIntFromTo(0, 5) <= 1) {
			// Randomly pick a point with collagen and move there
			int index = RandomHelper.nextIntFromTo(0, ecmWithCol.size() - 1);
			Object ecmAgent = ecmWithCol.get(index);
			GridPoint movePt = grid.getLocation(ecmAgent);
			grid.moveTo(this, movePt.getX(), movePt.getY());
		}
		else {
			int movement = 3;
			while (movement > 0) {
				// find neighbor with collagen or necrosis
				MooreQuery<Object> query = new MooreQuery(grid, this, 1, 1); // get neighboring agents with 1 grid unit
				Iterable<Object> iter = query.query(); // query the list of agents that are the neighbors
				List<Object> openNeighbor = new ArrayList<Object>();
				for (Object obj : iter){
					if (obj instanceof Necrosis){ // add to neighbor if necrotic
						openNeighbor.add(obj);
					}
					if (obj instanceof ECM && ((ECM) obj).getCollagen() > 0){ // add to neighbor if collagen > 0
						openNeighbor.add(obj);
					}
				}
				if (openNeighbor.size() == 0) { // if there's no necrotic or collagen neighbor
					// Randomly pick a point with collagen and move there
					int index = RandomHelper.nextIntFromTo(0, ecmWithCol.size() - 1);
					Object ecmAgent = ecmWithCol.get(index);
					GridPoint movePt = grid.getLocation(ecmAgent);
					grid.moveTo(this, movePt.getX(), movePt.getY());
				}
				else {
					int index = RandomHelper.nextIntFromTo(0, openNeighbor.size() - 1);
					Object randNeighbor = openNeighbor.get(index);
					GridPoint movePt = grid.getLocation(randNeighbor);
					grid.moveTo(this, movePt.getX(), movePt.getY());
					movement = movement - 1;
				}
			}
		}
	}
	
	// M1 BEHAVIOR
	public void m1Behave(Context<Object> context) {
		// omit M1 secretions
		
		// Find neutrophil and necrosis in this grid point
		MooreQuery<Object> query = new MooreQuery(grid, this, 1, 1); // get "cells" this macrophage is currently in
		Iterable<Object> iter = query.query(); // query the list of agents
		List<Object> necroticCell = new ArrayList<Object>();
		List<Object> neutrophilNgh = new ArrayList<Object>();
		for (Object obj : iter) {
			if (obj instanceof Necrosis) { // is current "cell" necrotic?
				necroticCell.add(obj);
			}
			else if (obj instanceof Neutrophil && ((Neutrophil) obj).getApoptosed() > 0) { // neutrophil in same "cell"?
				neutrophilNgh.add(obj);
			}
		}
		
		// Check on conditions to determine behavior
		if (linked == true) {
			// Macrophages will move towards its buddy. If there is none, its no longer linked			
			if (buddy != null) {
				GridPoint movePt = grid.getLocation(buddy);
				//System.out.println("Buddy: ( " + movePt.getX() + " , " + movePt.getY() + " )");
				grid.moveTo(this, movePt.getX(), movePt.getY());
			}
			else {
				linked = false;
			}
		} 
		else {
			if (necroticCell.size() > 0) {
				if (neutrophilNgh.size() > 0) {
					// choose randomly to each either necrosis (debris) or neutrophil
					if (RandomHelper.nextIntFromTo(0, 1) == 0) { // 0 = eat debris
						this.phagocytosis = phagocytosis - 1;
						// create secondary necrosis
						//Necrosis necrosis = new Necrosis(grid, 1, 0); // secondary = 1; age = 0
						//context.add(necrosis);
						//GridPoint gridpt = grid.getLocation(this);
						//grid.moveTo(necrosis, gridpt.getX(), gridpt.getY());
						// release ROS
						if (this.phagocytosis < -8) {
							if (RandomHelper.nextIntFromTo(0, 2) == 1) {
								macDie();
							}
						}
					} 
					else { // 1 = eat neutrophil
						this.phagocytosis = phagocytosis + 1;
						int index = RandomHelper.nextIntFromTo(0, neutrophilNgh.size() - 1);
						Object randNeutrophil = neutrophilNgh.get(index);
						context.remove(randNeutrophil);
					}
				}
				else { // eat necrosis (debris)
					this.phagocytosis = phagocytosis - 1;
					// create secondary necrosis
					//Necrosis necrosis = new Necrosis(grid, 1, 0); // secondary = 1; age = 0
					//context.add(necrosis);
					//GridPoint gridpt = grid.getLocation(this);
					//grid.moveTo(necrosis, gridpt.getX(), gridpt.getY());
					// release ROS
					if (this.phagocytosis < -8) {
						if (RandomHelper.nextIntFromTo(0, 2) == 1) {
							macDie();
						}
					}
				}
			}
			else { // if not necrosis at current location
				if (neutrophilNgh.size() > 0) {
					// eat neutrophil
					this.phagocytosis = phagocytosis + 1;
					int index = RandomHelper.nextIntFromTo(0, neutrophilNgh.size() - 1);
					Object randNeutrophil = neutrophilNgh.get(index);
					context.remove(randNeutrophil);
				}
				else {
					m1Migrate();
					if (RandomHelper.nextIntFromTo(0, 24) == 1) { // NetLogo orig: 5 * 4 + 5; chance M1 will polarize to M2
						stay = 1;
						phenotype = 1;
					}
				}
			}
		}
	}

	public void m1Migrate() {
		Context context = ContextUtils.getContext(this);
		List<Object> ecmEdgeAgent = ECM.getECMEdgeElems(context); // get list of ECM edge agents
		MooreQuery<Object> query = new MooreQuery(grid, this, 0, 0); // get "cells" this macrophage is currently in
		Iterable<Object> iter = query.query(); // query the list of agents
		List<Object> boundary = new ArrayList<Object>();
		for (Object obj : iter) {
			if (obj.equals(ecmEdgeAgent)) { // if queried agent is an ECM edge, add to boundary "list"
				boundary.add(obj);
			}
		}
		if (boundary.size() == 0) {
			if (growthFactors[14] > 0) { // check MCP (growth factor [14]) in current "patch" - currently unable due to
											// ODE - no info about spatial distribution
				// move up gradient 3x
				int moved = 0;
				while (moved < 4) {
					// detect MCP gradient and move up, increased moved
					GridPoint pt = grid.getLocation(this);
//						GridCellNgh<MCP> nghCreator = new GridCellNgh<MCP>(grid, pt, MCP.class, 2, 2);
//						List<GridCell<MCP>> gridCells = nghCreator.getNeighborhood(true);
//						GridPoint pointWithMostMCP = null;
//						int maxCount = -1;
//						for (GridCell<MCP> cell : gridCells) {
//							if (cell.size() > maxCount) {
//								pointWithMostMCP = cell.getPoint();
//								maxCount = cell.size();
//							}
//						}
//						double dX = pointWithMostMCP.getX() - pt.getX();								// calc x displacement
//						double dY = pointWithMostMCP.getY() - pt.getY();								// calc y displacement
//						double angle = SpatialMath.angleFromDisplacement(dX, dY);						// calc angle
//						grid.moveByVector(this, 1, angle);												// move 1 unit in dirxn of MCP gradient
					moved = moved + 1;

					// find neighboring apoptosed neutrophils
					pt = grid.getLocation(this); // update this Mac's location
					MooreQuery<Object> query2 = new MooreQuery(grid, this, 2, 2);					// find neighboring agents within 1.5 units
					Iterable<Object> iter2 = query2.query(); 										// query the list of agents
					List<Object> apopNeutrophil = new ArrayList<Object>();				
					for (Object obj : iter2) {
						if (obj instanceof Neutrophil && ((Neutrophil) obj).getApoptosed() > 0) {	// if queried agent is an apop. neu., add to list
							apopNeutrophil.add(obj);
						}
					}
					if (apopNeutrophil.size() > 0) {												// if there is an apop. neutrophil
						int index = RandomHelper.nextIntFromTo(0, apopNeutrophil.size() - 1);			// select a random one
						Object neutrophil = apopNeutrophil.get(index);
						GridPoint targetPt = grid.getLocation(neutrophil);							// get the neutrophil's location
						double dX = targetPt.getX() - pt.getX();									// calc x displacement
						double dY = targetPt.getY() - pt.getY();									// calc y displacement
						double disp = Math.hypot(dX, dY);											// calc total displacement
						dX = dX / disp;																// create unit vectors
						dY = dY / disp;
						grid.moveByDisplacement(this, (int) dX, (int) dY);
						moved = 5;
					}

					// find neighboring SSCs
					List<Object> buddyList = linkSSCbuddy();
					if (buddyList.size() > 0) {
						this.linked = true;
						this.buddy = buddyList.get(0);
					}
					else {
						this.linked = false;
						this.buddy = null;
					}
				}
			} 
			else {
				// move to "cell" with collagen greater than 0
				List<Object> ecms = ECM.getECM(context); // get all ECM agents
				List<Object> ecmWithCol = new ArrayList<Object>(); // create new array
				for (Object obj : ecms) { // search for ECM with collagen > 0
					if (((ECM) obj).getCollagen() > 0) {
						ecmWithCol.add(obj);
					}
				}
				int index = RandomHelper.nextIntFromTo(0, ecmWithCol.size() - 1);
				Object target = ecmWithCol.get(index);
				GridPoint movept = grid.getLocation(target);
				grid.moveTo(this, movept.getX(), movept.getY());
			}
		} 
		else {
			// move to "cell" with collagen greater than 0
			List<Object> ecms = ECM.getECM(context); // get all ECM agents
			List<Object> ecmWithCol = new ArrayList<Object>(); // create new array
			for (Object obj : ecms) { // search for ECM with collagen > 0
				if (((ECM) obj).getCollagen() > 0) {
					ecmWithCol.add(obj);
				}
			}
			int index = RandomHelper.nextIntFromTo(0, ecmWithCol.size() - 1);
			Object target = ecmWithCol.get(index);
			GridPoint movept = grid.getLocation(target);
			grid.moveTo(this, movept.getX(), movept.getY());
		}
	}
	
	// M1 EATS APOPTOTIC NEUTROPHILS
	public void m1ApopEat(Context<Object> context) {
		//Context context = ContextUtils.getContext(this);
		// secrete cytokines
		if (stay > 0) {
			stay = stay + 1;
			age = age - 1;
			if (stay > RandomHelper.nextIntFromTo(0, 2)) { // originally, NetLogo : random (2 * 1.5)
				stay = 0;
				this.phenotype = 3;
			}
		}
		else {
			MooreQuery<Object> query = new MooreQuery(grid, this, 0, 0); 					// get "cells" this macrophage is currently in
			Iterable<Object> iter = query.query(); 											// query the list of agents
			List<Object> apopNeutrophil = new ArrayList<Object>();				
			for (Object obj : iter) {
				if (obj instanceof Neutrophil && ((Neutrophil) obj).getApoptosed() > 0) {		// if queried agent is an apop. neu., add to list
					apopNeutrophil.add(obj);
				}
			}
			if (apopNeutrophil.size() > 0) {
				phagocytosis = phagocytosis + 1;
				int index = RandomHelper.nextIntFromTo(0, apopNeutrophil.size() - 1);
				Object neutrophil = apopNeutrophil.get(index);
				context.remove(neutrophil);
			}
			else {
				m1Migrate();
			}
			if (age > (1 + 8 * 4.3)) {
				stay = 1;
				age = 0;
			}
		}
	}
	
	//M1 EATS DEBRIS
	public void m1DebrisEat(Context<Object> context) {
		//Context context = ContextUtils.getContext(this);
		//secrete cytokine
		MooreQuery<Object> query = new MooreQuery(grid, this, 0, 0); // get "cells" this macrophage is currently in
		Iterable<Object> iter = query.query(); // query the list of agents
		List<Object> necroticCell = new ArrayList<Object>();
		for (Object obj : iter) {
			if (obj instanceof Necrosis) { // is current "cell" necrotic?
				necroticCell.add(obj);
			}
		}
		if (necroticCell.size() > 0) {
			phagocytosis = phagocytosis - 1;
			for (int i = 0; i < 2; i++) {
				//Necrosis necrosis = new Necrosis(grid, 1, 0);
				//context.add(necrosis);
				//GridPoint pt = grid.getLocation(this);
				//grid.moveTo(necrosis, pt.getX(), pt.getY());
			}
			// increase IL1 in "patch"
			// increase ROS within radius of 1.5
		}
		else {
			MooreQuery<Object> query2 = new MooreQuery(grid, this, 1, 1);
			Iterable<Object> iter2 = query2.query();
			List<Object> necroticNgh = new ArrayList<Object>();
			for (Object obj : iter2) {
				if (obj instanceof Necrosis) {
					necroticNgh.add(obj);
				}
			}
			if (necroticNgh.size() > 0) {
				int index = RandomHelper.nextIntFromTo(0, necroticNgh.size() - 1);
				Object necNeighbor = necroticNgh.get(index);
				GridPoint movePt = grid.getLocation(necNeighbor);
				grid.moveTo(this, movePt.getX(), movePt.getY());
			}
			else {
				nMigrate();
			}
		}
	}
	
	//M2 BEHAVIOR
	public void m2Behave(Context<Object> context) {
		// secrete cytokine
		int random =  RandomHelper.nextIntFromTo(0, 5);
		if (random == 1) {
			//macDie();
		}
	}
	
	// MACROPHAGE PROLIFERATION
	public void macProliferate() {
		Context context = ContextUtils.getContext(this);
		GridPoint gridpt = grid.getLocation(this);
		Macrophage macrophage = new Macrophage(space, grid, this.phenotype, 5, this.phagocytosis, this.linked, this.buddy); // set phenotype = parent, age = 5; phagocytosed = same as parent; linked = same as parent; buddy = same buddy as parent
		context.add(macrophage);
		//Object macrophage = this.clone(); // create copy of this agent (NetLogo hatch)
		grid.moveTo(macrophage, gridpt.getX(), gridpt.getY());		
	}
	
	// MACROPHAGE DEATH
	public void macDie() {
		Context context = ContextUtils.getContext(this);
		if (this.buddy != null) {
			//let the satellite buddy know this macrophage is gone - how?
		}
		context.remove(this);
	}
	
	// Link SSC buddy
	public List<Object> linkSSCbuddy() {
		MooreQuery<Object> query = new MooreQuery(grid, this, 1, 1); // get "cells" this macrophage is currently in
		Iterable<Object> iter = query.query(); // query the list of agents
		List<Object> sscList = new ArrayList<Object>();
		for (Object obj : iter) {
			if (obj instanceof SSC) {
				sscList.add(obj);
			}
		}
		List<Object> buddyList = new ArrayList<Object>();
		if (sscList.size() > 0) {
			int index = RandomHelper.nextIntFromTo(0, sscList.size() - 1);
			Object buddy = sscList.get(index);
			
			buddyList.add(buddy);
		}
		return buddyList;
	}
	
	public void nMigrate() {			// taken from Neutrophil.java
		Context context = ContextUtils.getContext(this);
		// At each step neutrophils should sense surroundings and move toward necrosis/high collagen
		MooreQuery<Object> query = new MooreQuery(grid, this, 2, 2); // get neighbors in a wider area
		Iterable<Object> iter = query.query(); // query the list of agents that are the neighbors
		List<Object> necrosisNeighbor = new ArrayList<Object>();
		List<Object> openNeighbor = new ArrayList<Object>();
		List<Object> highCollNeighbor = new ArrayList<Object>();
		for (Object neighbor : iter){
			if (neighbor instanceof Necrosis || neighbor instanceof ECM){
				openNeighbor.add(neighbor); // if any neighbors are necrotic or ecm, add to the list
			}
			if (neighbor instanceof Necrosis){ // necrotic neighbors
				necrosisNeighbor.add(neighbor);
			}
			if (neighbor instanceof ECM && ((ECM) neighbor).getCollagen() > 1){ // neighbors with collagen
				highCollNeighbor.add(neighbor);
			}
		}
		if (necrosisNeighbor.size() > 0 ){ // if there is necrosis move there
			int index = RandomHelper.nextIntFromTo(0, necrosisNeighbor.size()-1);
			Object randomNeighbor = necrosisNeighbor.get(index);
			GridPoint pt = grid.getLocation(randomNeighbor);
			grid.moveTo(this, pt.getX(), pt.getY());
		}
		else if (highCollNeighbor.size() > 0 ){ // if there is any collagen
			int index = RandomHelper.nextIntFromTo(0, highCollNeighbor.size()-1);
			Object randomNeighbor = highCollNeighbor.get(index);
			GridPoint pt = grid.getLocation(randomNeighbor);
			grid.moveTo(this, pt.getX(), pt.getY());	
		}
		else { // Otherwise just pick a random direction of ecm go to it
			int index = RandomHelper.nextIntFromTo(0, openNeighbor.size()-1);
			Object randomNeighbor = openNeighbor.get(index);
			GridPoint pt = grid.getLocation(randomNeighbor);
			grid.moveTo(this, pt.getX(), pt.getY());
		}				
		
	}
	
	// RETURN M1 AGENTS
	public static List<Object> getM1(Context<Object> context){ // Get a list of all the macrophages
		List<Object> m1Mac = new ArrayList<Object>(); // create a list of all the macrophage agents
		for (Object obj : context){ 
			if (obj instanceof Macrophage && ((Macrophage) obj).phenotype == 1){
				m1Mac.add(obj);
			}
		} 
		return m1Mac;
	}
		
	// RETURN M1 COUNT
	public int getM1Count() {		
		return m1num;
	}
	
	// RETURN M2 AGENTS
	public static List<Object> getM2(Context<Object> context){ // Get a list of all the macrophages
		List<Object> m2Mac = new ArrayList<Object>(); // create a list of all the macrophage agents
		for (Object obj : context){ 
			if (obj instanceof Macrophage && ((Macrophage) obj).phenotype == 2){
				m2Mac.add(obj);
			}
		} 
		return m2Mac;
	}
	
	// RETURN M2 COUNT
	public int getM2Count() {
		return m2num;
	}
	
	// RETURN RESIDENT MAC AGENTS
	public static List<Object> getMres(Context<Object> context){ // Get a list of all the macrophages
		List<Object> mresMac = new ArrayList<Object>(); // create a list of all the macrophage agents
		for (Object obj : context){ 
			if (obj instanceof Macrophage && ((Macrophage) obj).phenotype == 3){
				mresMac.add(obj);
			}
		} 
		return mresMac;
	}
	
	// RETURN RES MAC COUNT
	public int getMresCount() {
		return resMnum;
	}
	
	// RETURN M0 AGENTS (in the event that a phenotype was not assigned)
	public static List<Object> getM0(Context<Object> context){ // Get a list of all the macrophages
		List<Object> m0Mac = new ArrayList<Object>(); // create a list of all the macrophage agents
		for (Object obj : context){ 
			if (obj instanceof Macrophage && ((Macrophage) obj).phenotype == 0){
				m0Mac.add(obj);
			}
		} 
		return m0Mac;
	}
	
	// RETURN MACROPHAGE PHENOTYPE
	public int getPhenotype(){
		return phenotype;
	}
	
	// RETURN APOP EATING MACROPHAGES
	public static List<Object> getM1ae(Context<Object> context){ // Get a list of all the macrophages
		List<Object> aeMac = new ArrayList<Object>(); // create a list of all the macrophage agents
		for (Object obj : context){ 
			if (obj instanceof Macrophage && ((Macrophage) obj).phagocytosis > 0){
				aeMac.add(obj);
			}
		} 
		return aeMac;
	}
	
	public double getM1aeCount() {
		Context context = ContextUtils.getContext(this);
		double m1ae = getM1ae(context).size();
		return m1ae;
	}
	
	// RETURN DEBRIS EATING MACROPHAGES
	public static List<Object> getM1de(Context<Object> context){ // Get a list of all the macrophages
		List<Object> deMac = new ArrayList<Object>(); // create a list of all the macrophage agents
		for (Object obj : context){ 
			if (obj instanceof Macrophage && ((Macrophage) obj).phagocytosis < 0){
				deMac.add(obj);
			}
		} 
		return deMac;
	}
	
	public double getM1deCount() {
		Context context = ContextUtils.getContext(this);
		double m1de = getM1de(context).size();
		return m1de;
	}
	
}
