/**
 * 
 */
package muscle;

import java.util.ArrayList;
import java.util.List;

import repast.simphony.context.Context;
import repast.simphony.engine.schedule.ScheduledMethod;
import repast.simphony.query.space.grid.MooreQuery;
import repast.simphony.random.RandomHelper;
import repast.simphony.space.continuous.ContinuousSpace;
import repast.simphony.space.grid.Grid;
import repast.simphony.space.grid.GridPoint;
import repast.simphony.util.ContextUtils;

/**
 * @author kbc7vn
 *
 */
public class Neutrophil {

	// Kyle's code: neutrophils-own [phagocytosis speed stay taxisx taxisy apoptosed
	// age]

	// Neutrophil parameters
	private static Grid<Object> grid;
	private static ContinuousSpace<Object> space;
	public int phagocytosis; // neutrophil phagocytosis: starts at 0, 1+ converting to apoptotic neutrophil
								// after enough phagocytosis
	public int age; // neutrophil recruitment time- counter
	public int apoptosed; // neutrophil apoptosis counter

	/*
	 * // get growth factor numbers private static double[] growthFactors =
	 * GrowthFactors.getGrowthFactors(); // Recruiting growth factors private static
	 * double mcp = growthFactors[14]; private static double il1 = growthFactors[7];
	 * private static double il8 = growthFactors[8]; private static double ccl4 =
	 * growthFactors[12]; private static double cxcl1 = growthFactors[10]; private
	 * static double cxcl2 = growthFactors[9]; private static double gcsf =
	 * growthFactors[20]; // Deterring growth factors private static double lipoxins
	 * = growthFactors[22]; private static double mmp12 = growthFactors[19]; private
	 * static double resolvins = growthFactors[23]; private static double pge2 =
	 * growthFactors[27]; private static double lactoferrins = growthFactors[16];
	 * private static double il6 = growthFactors[13]; private static double il10 =
	 * growthFactors[21];
	 */

	public double numApopNeutrophils;
	public double numNeutrophils;
	public double totalNeutrophilCount;

	// Recruitment-promoting: mcp1 (14), il1 (7), il8 (8), ccl4 (12), cxcl1 (10),
	// cxcl2 (9)
	// Recruitment-deterring: lipoxins (22), mmp12 (19), resolvins (23), pge2 (27),
	// lactoferrins (16), il6 (13), il10 (21)

	public Neutrophil(Grid<Object> grid, ContinuousSpace<Object> space, int age, int phagocytosis, int apoptosed) {
		this.grid = grid;
		this.space = space;
		this.age = age;
		this.phagocytosis = phagocytosis;
		this.apoptosed = apoptosed;
	}
	
	
	@ScheduledMethod(start = 1, interval = 1, priority = 1, pick = 1)
	public void neutrophil_recruit() {// , Grid<Object> grid) { // recruitment/deter seesaw
		Context context = ContextUtils.getContext(this);
		double[] growthFactors = GrowthFactors.getGrowthFactors();
		// Recruiting growth factors
		double mcp = growthFactors[14];
		double il1 = growthFactors[7];
		double il8 = growthFactors[8];
		double ccl4 = growthFactors[12];
		double cxcl1 = growthFactors[10];
		double cxcl2 = growthFactors[9];
		double gcsf = growthFactors[20];
		// Deterring growth factors
		double lipoxins = growthFactors[22];
		double mmp12 = growthFactors[19];
		double resolvins = growthFactors[23];
		double pge2 = growthFactors[27];
		double lactoferrins = growthFactors[16];
		double il6 = growthFactors[13];
		double il10 = growthFactors[21];

		// determine how strong the recruitment force is
		double recruit = il1 + gcsf + mcp + ccl4 + cxcl2 + cxcl1 + il8;
		double deter = lipoxins + resolvins + mmp12 + lactoferrins + pge2 + il10 + il6 / 2;
		double differential = recruit - deter;

		List<Object> necrosis = Necrosis.getNecrosis(context);
		if (differential > 0) {
			int max = Fiber.origFiberNumber * 2;
			int neutro_prob = RandomHelper.nextIntFromTo(0, (int) ((differential / 2) + (differential / 4) - 1)/10);
			if (neutro_prob > 2) {
				neutro_prob = (2) + RandomHelper.nextIntFromTo(0, 1);
			}
			//System.out.print("neutro_prob: ");
			//System.out.print(neutro_prob);
			//System.out.print(" ");
			while (neutro_prob > 0 && getNeutrophilCount() < max) { // make some neutrophils
				// create_neutrophil();
				// find areas of necrosis to place neutrophil

				if (necrosis.size() != 0) {
					int index = RandomHelper.nextIntFromTo(0, necrosis.size() - 1);
					Object randomNecrosis = necrosis.get(index); // get the necrosis based on the random number chosen index
					
/*					if (getNeutrophilCount() <= 1) {
						//System.out.println("neutro");
						int i = 0;
						while (i<5) {
							GridPoint pt = grid.getLocation(randomNecrosis); // Get the necrosis location
							Neutrophil newNeutrophil = new Neutrophil(grid, space, 0, 0, 0); // change to neutrophil
							context.add((Neutrophil) newNeutrophil);
							grid.moveTo(newNeutrophil, pt.getX(), pt.getY());
							i = i + 1;
							//System.out.println(getNeutrophilCount());
						}
					
					} */
					
					try { // within the list
						//System.out.print(necrosis.size());
						//System.out.print(" ");
						//System.out.print(index);
						//System.out.print(" ");

						GridPoint pt = grid.getLocation(randomNecrosis); // Get the necrosis location
						Neutrophil newNeutrophil = new Neutrophil(grid, space, 0, 0, 0); // change to neutrophil
						context.add((Neutrophil) newNeutrophil);
						grid.moveTo(newNeutrophil, pt.getX(), pt.getY());
						//System.out.print(pt);

					} catch (Exception e) {
						e.printStackTrace();
					}

				}
				neutro_prob = neutro_prob - 1;
			}
		}
	}

	@ScheduledMethod(start = 2, interval = 1, priority = 2)
	public void step() {
		Context context = ContextUtils.getContext(this);

		// Age neutrophil
		this.age = age + 1;

		if (this.apoptosed == 0) {
			// Find areas with necrosis

			// Find necrosis in this grid point
			MooreQuery<Object> query = new MooreQuery<Object>(grid, this, 0, 0); // get the location this neutrophil is
			// currently in
			Iterable<Object> iter = query.query(); // query the list of agents
			List<Object> necroticCell = new ArrayList<Object>();
			for (Object obj : iter) {
				if (obj instanceof Necrosis) { // is current "cell" necrotic?
					necroticCell.add(obj);
				}
			}

			if (necroticCell.size() > 0) {
				int index = RandomHelper.nextIntFromTo(0, necroticCell.size() - 1);
				Object necrosisRemove = necroticCell.get(index);
				GridPoint pt = grid.getLocation(necrosisRemove);
				// neutrophil_phagocytose(context, pt);
				this.phagocytosis = phagocytosis + 1;

				// Include when growth factor secretions are added:
				/*
				 * ifelse NO < 8 [ ask patches in-radius 2 [set ROS ROS + .5] ] [ ask patches
				 * in-radius 1 [set ROS ROS + .5] ]
				 * 
				 * ask patch-here [set necrotic necrotic + 1 // Note from Katherine - I am
				 * confused about this line, why would necrosis increase here? 
				 * set TNF TNF + 1
				 * ;They definitely secrete these three factors, just unsure if it happens in
				 * sterile inflammation, and to what extent 
				 * set IL1 IL1 + 1 set IL6 IL6 + 1 
				 * set IL8 IL8 + 1 set MCP MCP + 1.1 * 8.9 ] 
				 * ;They also make these cyto and chemokines
				 * set CCL4 CCL4 + 1 
				 * set CCL3 CCL3 + 1 
				 * set CXCL1 CXCL1 + 1 set CXCL2 CXCL2 + 1 set IFN IFN + 1 * 3.5
				 */
			} else {
				// neutrophil_migrate(context);
			}

			int index = RandomHelper.nextIntFromTo(0, 3 + 3 - 1); // subtract 1 to account for Netlogo to Repast
			if (index <= phagocytosis) { // After enough phagocytoses, it either becomes apoptotic or leaves
				this.apoptosed = apoptosed + 1;
				this.age = 0;
				// there is a 1/3 chance the neutrophil will now die
				int chance_die = RandomHelper.nextIntFromTo(0, 2);
				if (chance_die == 1) {
					context.remove(this);
				}
			}
		} else { // This is when apoptosed = 1, the promote M1 macrophages and deter neutrophils
			/*
			 * Include when growth factor secretions are added: set Lactoferins Lactoferins
			 * + 5 set HGF HGF + 1 set VEGF VEGF + 1 if CCL3 > 0 [set CCL3 CCL3 - 1]
			 */
			int chance_die = RandomHelper.nextIntFromTo(0, 2);
			if (chance_die == 1) {
				context.remove(this);
			}
			
		}

		// Once a neutrophil is 8 hours old, it has a chance of leaving, more likely the
		// longer its around
		if (this.age > 3 + 5 * 2) {
			int chance_die = RandomHelper.nextIntFromTo(0, 4);
			if (chance_die == 1) {
				context.remove(this);
			}
		}

	}

	// ADD NEW NEUTROPHILS TO NECROTIC TISSUE
	public void create_neutrophil() {
		Context context = ContextUtils.getContext(this);

		// find areas of necrosis to place neutrophil
		List<Object> necrosis = Necrosis.getNecrosis(context);
		if (necrosis.size() != 0) {
			int index = RandomHelper.nextIntFromTo(0, necrosis.size() - 1); // instead of random number, should they
																			// appear at edge necrosis first?
			Object randomNecrosis = necrosis.get(index); // get the necrosis based on the random number chosen index
															// within the list
			GridPoint pt = grid.getLocation(randomNecrosis); // Get the necrosis location

			Neutrophil newNeutrophil = new Neutrophil(grid, space, 0, 0, 0); // change to neutrophil
			context.add(newNeutrophil);
			grid.moveTo(newNeutrophil, pt.getX(), pt.getY());
		}

		/*
		 * Include once we set up growth factor secretion growth factor secretion from
		 * neutrophil creation set Lactoferins Lactoferins + 2 ;neutrophils release
		 * Lactoferins when they exit the blood set Azurocidin Azurocidin + 2 set LL37
		 * LL37 + 2 set CathepsinG CathepsinG + 2
		 */

	}

	public void neutrophil_migrate() {
		Context  context = ContextUtils.getContext(this);
		// At each step neutrophils should sense surroundings and move toward
		// necrosis/high collagen
		MooreQuery<Object> query = new MooreQuery(grid, this, 2, 2); // get neighbors in a wider area
		Iterable<Object> iter = query.query(); // query the list of agents that are the neighbors
		List<Object> necrosisNeighbor = new ArrayList<Object>();
		List<Object> openNeighbor = new ArrayList<Object>();
		List<Object> highCollNeighbor = new ArrayList<Object>();
		List<Object> ECMNeighbor = new ArrayList<Object>();
		for (Object neighbor : iter) {
			if (neighbor instanceof Necrosis) {
				necrosisNeighbor.add(neighbor); // if any neighbors are necrotic, add to the list
			} else if (neighbor instanceof ECM && ((ECM) neighbor).getCollagen() > 1) { // neighbors with collagen
				highCollNeighbor.add(neighbor);
			} else if (neighbor instanceof ECM) { // neighbors with collagen
				ECMNeighbor.add(neighbor);
			} else {
				openNeighbor.add(neighbor);
			}
		}
		// If there is necrosis move there
		if (necrosisNeighbor.size() > 0) {
			int index = RandomHelper.nextIntFromTo(0, necrosisNeighbor.size() - 1);
			Object randomNeighbor = necrosisNeighbor.get(index);
			GridPoint pt = grid.getLocation(randomNeighbor);
			grid.moveTo(this, pt.getX(), pt.getY());
		}
		// If there is no necrosis, move towards an area with collagen
		else if (highCollNeighbor.size() > 0) {
			int index = RandomHelper.nextIntFromTo(0, highCollNeighbor.size() - 1);
			Object randomNeighbor = highCollNeighbor.get(index);
			GridPoint pt = grid.getLocation(randomNeighbor);
			grid.moveTo(this, pt.getX(), pt.getY());
		}
		// Otherwise just pick a random direction to move
		else {
			int index = RandomHelper.nextIntFromTo(0, openNeighbor.size() - 1);
			Object randomNeighbor = openNeighbor.get(index);
			GridPoint pt = grid.getLocation(randomNeighbor);
			grid.moveTo(this, pt.getX(), pt.getY());
		}

	}

	public void neutrophil_phagocytose(Context<Object> context, GridPoint pt) {
		this.phagocytosis = phagocytosis + 1;
		// Remove necrotic element
		context.remove(pt);
		// Replace with ECM element
		ECM newECM = new ECM(grid, 0.1, 0, 0);
		context.add(newECM);
		// Move ECM to where the necrosis was removed
		grid.moveTo(newECM, pt.getX(), pt.getY());
	}

	// Get the total apoptotic neutrophil count
	public static List<Object> getApoptosed(Context<Object> context) {
		//Context context = ContextUtils.getContext(this);
		List<Object> apopNeutrophils = new ArrayList<Object>();
		for (Object obj : context) {
			if (obj instanceof Neutrophil && ((Neutrophil) obj).apoptosed >= 1) {
				apopNeutrophils.add(obj);
			}
		}
		return apopNeutrophils;
	}
	
	public double getApoptosedCount() {
		Context context = ContextUtils.getContext(this);
		numApopNeutrophils = getApoptosed(context).size();
		return numApopNeutrophils;
	}

	// Get the total (non-apoptotic) neutrophil count
	public static List<Object> getNeutrophils(Context<Object> context) {
		//Context context = ContextUtils.getContext(this);
		List<Object> neutrophils = new ArrayList<Object>();
		for (Object obj : context) {
			if (obj instanceof Neutrophil && ((Neutrophil) obj).apoptosed == 0) {
				neutrophils.add(obj);
			}
		}
		return neutrophils;
	}
	
	public double getNeutrophilCount() {
		Context context = ContextUtils.getContext(this);
		numNeutrophils = getNeutrophils(context).size();
		return numNeutrophils;
	}
	
	
	public double getTotalNeutrophilCount() {
		totalNeutrophilCount = getNeutrophilCount() + getApoptosedCount();
		return totalNeutrophilCount;
	}

	public double getApoptosed() {
		return apoptosed;
	}

	public double getAge() {
		return age;
	}

	public double getPhagocytosis() {
		return phagocytosis;
	}

}
