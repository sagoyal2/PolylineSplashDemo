/** 
 * 	
 * 
 * 	PolylineSplashDemo
 * 	@author Samaksh (Avi) Goyal, 6/27/2021
 * 
 * 	Core properties:
 * 	Constant Volume Constraint
 * 	Converservation of Momentum Constraint
 * 
 * 	Syringe like Inflation / Deflation
 * 
 * 	Central Update:
 * 	MCF
 * 	project onto volume + mom
 *  Laplacian based Reprarameterization
 *
 * 
 * 
 * 
 * 
 */

float BRUSH_RADIUS = 50.0;
float MESH_RESOLUTION = 50;
float INITIAL_DROPLET_RADIUS = 100;

boolean NORMAL_FLAG = false;
boolean CURVATURE_FLAG = false;
boolean WEIGHT_FLAG = false;
boolean WEIGHT_MODE = false;
boolean MCF_FLAG = false;
int mcf_clock = 100;
boolean GEODESIC_FLAG = false;
boolean IMAGE_FLAG = false;
boolean DRIG_FLAG = false;
boolean VOLUME_FLAG = false;
boolean DEPTH_FLAG = false;
boolean	SHOW_FUTURE = false;
boolean WITH_MOM_CONSTRAINT = false;

float initial_area = -1.0;
SplashBrush brush;
SplashBrush brish;
PolylineSplash my_splash;
LinkedList<PolylineSplash> undo_splash;
ArrayList<SplashBrush> pins  = new ArrayList<SplashBrush>();
ArrayList<SplashBrush> rigs  = new ArrayList<SplashBrush>();
SplashBrush current_rig;

// Drivers
/////////////////////////////////////////////////////////////////
void setup(){
	size(1280, 1024);  
	smooth(8);
	reset();
}

void reset(){ 

	my_splash 		= new PolylineSplash(width, height, MESH_RESOLUTION, INITIAL_DROPLET_RADIUS);
	undo_splash 	= new LinkedList<PolylineSplash>();

	initial_area = my_splash.getArea();
	println("initial_area: " + initial_area);
	pins.clear();
}


void draw(){
	background(255);


	fill(0);
	StringBuilder controls = new StringBuilder();
	controls.append("Radius[mouse wheel +/-], ");
	controls.append("undo[z], ");
	controls.append("show_normals[n], ");
	controls.append("label_curvature[k], ");
	controls.append("weight[w], ");
	controls.append("adjust weight[a], ");
	controls.append("MCF[m], ");
	controls.append("geodesic[g], ");
	controls.append("image[i], ");
	controls.append("drawrig[d], ");
	controls.append("depth[u], ");
	controls.append("add pin[p], ");
	controls.append("clear pin[c], ");
	String all_controls = controls.toString();
	text(all_controls, 5, 10);
	text("#UNDOS: " + undo_splash.size(), 5, 25);

	// Visualize
	drawBrushes();
	drawRigs();
	my_splash.refineMesh();
	// my_splash.fixDepth(); - fix this!
	my_splash.viewPoints(SHOW_FUTURE, WITH_MOM_CONSTRAINT);

	if(NORMAL_FLAG){
		my_splash.drawPointNormals();
		fill(52, 125, 235);
		text("NORMAL_FLAG ON - press[n] to remove", 5, 40); 
	}
	if(CURVATURE_FLAG){
		my_splash.labelCurvature();
		fill(237, 164, 245);
		text("CURVATURE_FLAG ON - press[k] to remove", 5, 55);
	}
	if(WEIGHT_FLAG){
		my_splash.showWeight();
		fill(237, 159, 50);
		text("WEIGHT_FLAG ON - press[w] to remove", 5, 70);
	}
	if(WEIGHT_MODE){
		fill(237, 159, 50);
		text("WEIGHT_ADJUSTMENT_MODE ON - press[a] to remove", 5, 85);
	}
	if(MCF_FLAG){
		if(mcf_clock < 5){
			fill(0);
			text("Applying Mean Curvature Flow on iteration i: " + mcf_clock, 5, 100);
			my_splash.mcf(mouseX, mouseY, BRUSH_RADIUS, GEODESIC_FLAG);
			mcf_clock++;
		}
		else{
			saveSplashState();
			MCF_FLAG = false;
			mcf_clock = 0;
		}

		/**
		 * NOTE THIS ALLOWS FOR FINE GRAIN CONTROL - small changes in mean curvature 
		 * ie, if the filter is NOT present then we get global curvature fix,
		 * if filter IS present then we get small detailed changes
		 * **/
		if(!GEODESIC_FLAG){

			applyWaterEffects(false);


			// // Make mesh spacing uniform
			// my_splash.reParameterize();

			// // Move onto constraint
			// my_splash.projectVolumePositions(initial_area);
		}
	}
	if(GEODESIC_FLAG){
		//my_splash.showGeodesic(mouseX, mouseY, BRUSH_RADIUS);
		fill(46, 166, 50);
		text("GEODESIC_FLAG ON - press[g] to remove", 5, 115);
	}
	if(IMAGE_FLAG){
		fill(235, 64, 52);
		text("IMAGE_FLAG ON - press[i] to remove", 5, 130);

		// dashed line at center:
		for(int i = 0; i < width; i+=10){
			ellipse(i, height/2, 2, 2);
		}
	}
	if(DRIG_FLAG){
		fill(172, 167, 176);
		text("DRIG_FLAG ON - press[d] to remove", 5, 145);
	}
	 if(VOLUME_FLAG){
	 	//my_splash.showVolume();
		fill(48, 174, 217);
		text("VOLUME_FLAG ON - press[v] to remove", 5, 160);
	}
	if(DEPTH_FLAG){
		// my_splash.showDepth();
		fill(252, 132, 3);
		text("DEPTH_FLAG ON - press[u] to remove", 5, 175);
	}
	fill(0);
}

void keyPressed(){
	if (key == 'r'){ 
		reset();
	}
	if (key == 'n'){
		NORMAL_FLAG = !NORMAL_FLAG;
	}
	if (key == 'k'){
		CURVATURE_FLAG = !CURVATURE_FLAG;
	}
	if (key == 'm'){
		MCF_FLAG = true;
	}
	if (key == 'w'){
		WEIGHT_FLAG = !WEIGHT_FLAG;
	}
	if (key == 'a'){
		WEIGHT_MODE = !WEIGHT_MODE;
	}
	if (key == 'z'){//undo
		if(DRIG_FLAG && rigs.size() > 0){
			rigs.remove(rigs.size() - 1);
		}
		else if (undo_splash.size() > 0) {
			my_splash = undo_splash.pollFirst();
		}
	}
	if (key == 'g'){
		GEODESIC_FLAG = !GEODESIC_FLAG;
	}
	if (key == 'i'){
		IMAGE_FLAG = !IMAGE_FLAG;
	}
	if(key == 'd'){
		DRIG_FLAG = !DRIG_FLAG;

		if(DRIG_FLAG){
			//set up rigs (pins)
		}else{
			//deform to make rig
			boolean with_rigs = true;
			boolean show_future = false;
			boolean with_mom_constraint = WITH_MOM_CONSTRAINT;
			boolean final_position = true;
			solveAndDeform(with_rigs, show_future, with_mom_constraint, final_position);
			rigs.clear();
		}
	}
	if(key == ENTER){
		SHOW_FUTURE = !SHOW_FUTURE;
		if(WITH_MOM_CONSTRAINT){
			WITH_MOM_CONSTRAINT = false;
		}
		// if(DRIG_FLAG){
		// 	boolean with_rigs = true;
		// 	boolean show_future = true;
		// 	boolean with_mom_constraint = false;
		// 	solveAndDeform(with_rigs, show_future, with_mom_constraint);
		// }
	}
	if(key == TAB){
		WITH_MOM_CONSTRAINT = !WITH_MOM_CONSTRAINT;
		if(DRIG_FLAG && SHOW_FUTURE){
			boolean with_rigs = true;
			boolean show_future = true;
			boolean with_mom_constraint = true;
			boolean final_position = false;
			solveAndDeform(with_rigs, show_future, with_mom_constraint, final_position);	
		}
	}
	if(key == 'v'){
		VOLUME_FLAG = !VOLUME_FLAG;
	}
	if(key == 'u'){
		DEPTH_FLAG = !DEPTH_FLAG;
	}
	if (key == 'p') {// pin/unpin
		PVector pinPos = new PVector(mouseX, mouseY, 0);
		// KelvinBrush pinBrush = new KelvinBrush(eps, nu);
		SplashBrush pinBrush = new SplashBrush(mouseX, mouseY, BRUSH_RADIUS);
		pinBrush.setPosition(pinPos);
		pins.add(pinBrush);
	}
	if (key == 'c') {// clear closest pin
		removeClosestPin();
	}
}

void mouseWheel(MouseEvent event) {

	float e = event.getCount();
	float scale = (e > 0.0)? 1.1:0.9;

	if (WEIGHT_MODE){
		reweight(scale);
	}

	else{

		if(VOLUME_FLAG){
			//QUOKKA the 0.5, 
			// we need to figure out a way to scale the amount of volume we want to increase
			// we might also want to apply some kinda of geodesic filter on where we insert the volume
			// we can't directly place brush inside of the volume in 3D, something to think about
			float volume_direction_scale = (scale > 1)? 1.0: -1.0;
			volume_direction_scale *= 0.1*BRUSH_RADIUS;
			changeVolume(volume_direction_scale, GEODESIC_FLAG);

			applyWaterEffects(true);

			// //update new area
			// initial_area = my_splash.getArea();

			// // Move onto constraint
			// my_splash.projectVolumePositions(initial_area);
			// // Make mesh spacing uniform
			// my_splash.reParameterize();

		}else{
			BRUSH_RADIUS *= scale;
		}
	}

}


/////////////////////////////////////////////////////////////////
void drawBrushes() {
	drawBrush(mouseX, mouseY);
	if(IMAGE_FLAG){
		drawBrush(mouseX, height - mouseY);
	}

	noFill();
	stroke(0, 255, 0);
	for (SplashBrush pin : pins) {
		// ellipse((float)pin.position.x, (float)pin.position.y, pin.radius, pin.radius);
		circle(pin.position.x, pin.position.y, 2*pin.radius);
	}

	drawRigs();
}

void drawRigs(){
	for (SplashBrush rig : rigs){
		stroke(255, 128, 89);
		circle(rig.position.x, rig.position.y, 2*rig.radius);
		circle(rig.position.x + rig.displacement.x, rig.position.y + rig.displacement.y, 2*rig.radius);

		line(rig.position.x, rig.position.y, rig.position.x + rig.displacement.x, rig.position.y + rig.displacement.y);

		if(WITH_MOM_CONSTRAINT){
			stroke(52, 174, 235);
			line(rig.position.x, rig.position.y, rig.position.x + rig.mom_conv_disp.x, rig.position.y + rig.mom_conv_disp.y);
		}

	}
}

void drawBrush(float x, float y)
{
	noFill();
	if(WEIGHT_MODE){
		stroke(237, 159, 50);
	}
	else{
		stroke(0,0,0);
	}
	circle(x, y, 2*BRUSH_RADIUS);
}

void mousePressed() {
	if(DRIG_FLAG){
		current_rig = new SplashBrush(mouseX, mouseY, BRUSH_RADIUS);
		//rigs.add(current_rig);
	}
	else{
		brush = new SplashBrush(mouseX, mouseY, BRUSH_RADIUS);
		if(IMAGE_FLAG){
			brish = new SplashBrush(mouseX, height - mouseY, BRUSH_RADIUS);
		}
	}

	saveSplashState();
}

void saveSplashState() {
	PolylineSplash previous = new PolylineSplash(my_splash);
	undo_splash.addFirst(previous);

	if (undo_splash.size() > 100) undo_splash.removeLast();
}

void mouseDragged() {

	if(DRIG_FLAG){
		PVector new_position = new PVector(mouseX, mouseY, 0);
		current_rig.setDisplacementBasedOnNewPosition(new_position);		

		noFill();
		stroke(255, 128, 89);
		SplashBrush rig = current_rig;
		circle(rig.position.x, rig.position.y, 2*rig.radius);
		circle(rig.position.x + rig.displacement.x, rig.position.y + rig.displacement.y, 2*rig.radius);
		line(rig.position.x, rig.position.y, rig.position.x + rig.displacement.x, rig.position.y + rig.displacement.y);

		boolean with_rigs = true;
		boolean show_future = SHOW_FUTURE;
		boolean with_mom_constraint = WITH_MOM_CONSTRAINT;
		boolean final_position = false;
		solveAndDeform(with_rigs, show_future, with_mom_constraint, final_position);

		// do I need to do QUOKKA?
		// current_rig.setPosition(new_position);
	}
	else{
		PVector new_position = new PVector(mouseX, mouseY, 0);
		brush.setDisplacementBasedOnNewPosition(new_position);

		PVector new_position_2 =  new PVector(mouseX, height - mouseY, 0);
		if(IMAGE_FLAG){
			brish.setDisplacementBasedOnNewPosition(new_position_2);
		}

		boolean with_rigs = false;
		boolean show_future = false;
		boolean with_mom_constraint = false;
		boolean final_position = true;
		solveAndDeform(with_rigs, show_future, with_mom_constraint, final_position);

		// Slide brush to newP:
		brush.setPosition(new_position);

		if(IMAGE_FLAG){
			brish.setPosition(new_position_2);
		}
	}
}

void mouseReleased(){
	if(DRIG_FLAG){
		rigs.add(current_rig);
		getMass();
		current_rig.getMomentum();
	}
}

void removeClosestPin() 
{
	if (pins.isEmpty()) return;

	// FIND PIN CLOSEST TO MOUSE:
	PVector p = new PVector(mouseX, mouseY, 0);
	int   minIndex = -1;
	float minDist  = Float.MAX_VALUE;
	for (int k=0; k<pins.size(); k++) {
		PVector pk = pins.get(k).position;
		float   dk = p.dist(pk);//sloth sqrt
		if (dk < minDist) {
			minDist = dk;
			minIndex = k;
		}
	}

	// REMOVE minIndex PIN:
	pins.remove(minIndex);
}


void solveAndDeform(boolean with_rigs, boolean show_future, boolean with_mom_constraint, boolean final_position)
{
	try {
		FastPinConstraintSolver2D solver;
		PVector norm_total_momentum = new PVector();
		float total_mass_sq = 0.0;

		if(with_rigs){
			if(!final_position){
				solver = new FastPinConstraintSolver2D(current_rig);
			}else{
				solver = new FastPinConstraintSolver2D();
			}

			for (SplashBrush rig : rigs) solver.addRig(rig);

			//if we have momentum constraints we will precompute sum(mi*vi)/mTm
			if(with_mom_constraint){
				for(SplashBrush rig : rigs){ 
					norm_total_momentum.add(rig.momentum);
					total_mass_sq += sq(rig.mass);
				}
				norm_total_momentum.div(total_mass_sq);

				// println("norm_total_momentum.x: " + norm_total_momentum.x + " norm_total_momentum.y: " + norm_total_momentum.y);
			}
		}
		else{
			solver = new FastPinConstraintSolver2D(brush);
			for (SplashBrush pin : pins) solver.addPin(pin);
		}

		//Set and solve u = K f, where here we are solving for f, after possibly updating u to be momemtum conserving
		solver.solve(with_rigs, with_mom_constraint, norm_total_momentum);

		if (!(GEODESIC_FLAG || DEPTH_FLAG)) {
			if(!final_position){
				my_splash.startFutureSplash();
				for (PVector p0 : my_splash.future_splash)  solver.deform(p0);
			}else{
				for (PVector p0 : my_splash.splash)  solver.deform(p0);

				boolean update_volume = true;
				applyWaterEffects(update_volume);
			}
		} 
		else {
			//int index = my_splash.indexOfClosestPoint(brush.position);
			//my_splash.getGeodesic(index);
			// my_splash.getGeodesic(brush.position.x, brush.position.y, BRUSH_RADIUS);
			ArrayList<Float> g = null;
			float maximum_g = Float.MIN_VALUE;

			if(GEODESIC_FLAG){
				my_splash.getGeodesic(brush.position.x, brush.position.y, BRUSH_RADIUS);
				g = my_splash.geodesic;
			}
			else if(DEPTH_FLAG){
				my_splash.getDepth();
				g = my_splash.depth;
			}

			// scale geodesic/ depth weight
			for (int k = 0; k < g.size(); k++){
				if(g.get(k) > maximum_g){
					maximum_g = g.get(k);
				}
			}

			// Deform with geodesic filter:
			for (int k=0; k<my_splash.splash.size(); k++) {
				PVector pk = my_splash.splash.get(k);
				PVector u  = solver.getDeformDisplacement(pk);
				// float   r  = pk.dist(brush.position);
				// float   re = (float)Math.sqrt(r*r + BRUSH_RADIUS*BRUSH_RADIUS);
				// float   gr = (g.get(k))/re;
				// float   a  = 2; // ** geodesic scale parameter **
				// if (gr > a) {
				// 	u.mult(1/(1+(gr-a)*(gr-a)));
				// }
				if(GEODESIC_FLAG){
					u.mult(pow((maximum_g - g.get(k))/(maximum_g), 2));
				}
				if(DEPTH_FLAG){
					u.mult(pow((g.get(k))/(maximum_g), 2));
				}
				pk.add(u);
			}

			boolean update_volume = true;
			applyWaterEffects(update_volume);
		}

	}
	catch(Exception e) {// singular matrix
		System.out.println("Matrix is singular: reverting to simple brush.");
		pins.clear();
		// brush.deform(P);
	}
}


void getMass(){
	
	int total_hits = 0;
	int offset = 4;

	// loop over all pixel locations on screen
	for(int i = 1; i < width; i+=offset){
		for(int j = 1; j < height; j+=offset){
			PVector point = new PVector(i, j, 0.0);
			boolean hit = my_splash.getIntersections(point);

			if(hit){
				total_hits++;

				//intersection successfull call getA2D_3x3Block (K matrix) on point
				float[] aij = new float[9];
				current_rig.getA2D_3x3Block(point, aij);

				//calculate forbenius norm
				float f_norm = 0.0;
				for (int i3=0; i3<3; i3++) {
          for (int j3=0; j3<3; j3++) {
            f_norm += sq(aij[3*i3+j3]);
          }
        }
        f_norm = sqrt(f_norm);

        //update brush/rig mass
        current_rig.mass += offset*offset*f_norm;
			}
		}
	}

	current_rig.mass /= total_hits;

	println("rig has mass: " + current_rig.mass);
}

void reweight(float weigth_scale){
	for(int i = 0; i < my_splash.splash.size(); i++){
		PVector p = my_splash.splash.get(i);
		float sqr_dist = (p.x - mouseX)*(p.x - mouseX) + (p.y - mouseY)*(p.y - mouseY);
		if(sqr_dist < BRUSH_RADIUS*BRUSH_RADIUS){ 
			float prior_weight = my_splash.weight.get(i);

			if(prior_weight > 100.0){
				prior_weight = 100.0;
			}
			my_splash.weight.set(i, weigth_scale*prior_weight);
		}
	}
}

void changeVolume(float direction, boolean with_geodesic){
	saveSplashState();

	if(with_geodesic){
		my_splash.getGeodesic(mouseX, mouseY, BRUSH_RADIUS);
	}

	float mu = 0.9;
	float nu = 0.3;
	float a = 1.0/(4.0*PI*mu);
	float b = a/(4.0*(1.0-nu));
	float eps = 0.01;
	float scaling = 200;

	for(int i = 0; i < my_splash.splash.size(); i++){
		PVector p = my_splash.splash.get(i);
		PVector r = new PVector(p.x - mouseX, p.y - mouseY, 0.0);
		float r_e = sqrt(r.magSq() + pow(eps,2));

		float geo_scaling = 1.0;
		float coef = 2*(b - a)*(1.0/pow(r_e,2) + (pow(eps, 2))/(2*pow(r_e, 4)));
		
		if(with_geodesic){
			geo_scaling *= my_splash.geodesic.get(i);
		}

		p.add(PVector.mult(r, geo_scaling*coef*scaling*direction));
	}
}

void applyWaterEffects(boolean reinitialize_area){

	if(reinitialize_area){
		initial_area = my_splash.getArea();
	}

	// Make mesh uniformally spaced
	my_splash.reParameterize();

	// Move onto volume constraint
	my_splash.projectVolumePositions(initial_area);

	//mean curvature flow
	my_splash.mcf(mouseX, mouseY, BRUSH_RADIUS, GEODESIC_FLAG);
}




























