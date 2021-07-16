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

float initial_area = -1.0;
SplashBrush brush;
SplashBrush brish;
PolylineSplash my_splash;
LinkedList<PolylineSplash> undo_splash;
ArrayList<SplashBrush> pins  = new ArrayList<SplashBrush>();

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
	my_splash.refineMesh();
	// my_splash.fixDepth(); - fix this!
	my_splash.viewPoints(DRIG_FLAG);


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
		my_splash.showGeodesic(mouseX, mouseY, BRUSH_RADIUS);
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
		fill(48, 174, 217);
		text("VOLUME_FLAG ON - press[v] to remove", 5, 160);
	}
	if(DEPTH_FLAG){
		my_splash.showDepth();
		fill(252, 132, 3);
		text("DEPTH_FLAG ON - press[u] to remove", 5, 175);
		//DEPTH_FLAG = false;
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
		if (undo_splash.size() > 0) {
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
			my_splash.startFutureSplash();
		}else{
			my_splash.projectToFuture();

			int iteration = 3;
			while(iteration > 0){
				// // Make mesh spacing uniform
				// my_splash.reParameterize();

				// // Move onto constraint
				// my_splash.projectVolumePositions(initial_area);


				applyWaterEffects(false);

				iteration--;
			}

			my_splash.setWeightToNeutral();
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

	brush = new SplashBrush(mouseX, mouseY, BRUSH_RADIUS);
	
	if(IMAGE_FLAG){
		brish = new SplashBrush(mouseX, height - mouseY, BRUSH_RADIUS);
	}

	saveSplashState();
}

void saveSplashState() {
	PolylineSplash previous = new PolylineSplash(my_splash);
	undo_splash.addFirst(previous);
}

void mouseDragged() {

	PVector new_position = new PVector(mouseX, mouseY, 0);
	brush.setDisplacementBasedOnNewPosition(new_position);

	PVector new_position_2 =  new PVector(mouseX, height - mouseY, 0);
	if(IMAGE_FLAG){
		brish.setDisplacementBasedOnNewPosition(new_position_2);
	}

	// if(DRIG_FLAG){
	// 	drawFuture();
	// }
	// else{
	// 	// move points within brush somehow
	// 	deform();
	// }

	solveAndDeform();

	// Slide brush to newP:
	brush.setPosition(new_position);

	if(IMAGE_FLAG){
		brish.setPosition(new_position_2);
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


void solveAndDeform()
{
	try {
		FastPinConstraintSolver2D solver = new FastPinConstraintSolver2D(brush);
		for (SplashBrush pin : pins) solver.addPin(pin);
		solver.solve();

		//if (!GEODESIC_FLAG) {
			for (PVector p0 : my_splash.splash)  solver.deform(p0);
		//} 
		// else {
		// 	int index = my_splash.indexOfClosestPoint(brush.position);
		// 	my_splash.getGeodesic(index);
		// 	ArrayList<Float> g = my_splash.geodesic;
		// 	// Deform with geodesic filter:
		// 	for (int k=0; k<my_splash.splash.size(); k++) {
		// 		PVector pk = my_splash.splash.get(k);
		// 		PVector u  = solver.getDeformDisplacement(pk);
		// 		float   r  = pk.dist(brush.position);
		// 		float   re = (float)Math.sqrt(r*r + BRUSH_RADIUS*BRUSH_RADIUS);
		// 		float   gr = (g.get(k))/re;
		// 		float   a  = 2; // ** geodesic scale parameter **
		// 		if (gr > a) {
		// 			u.mult(1/(1+(gr-a)*(gr-a)));
		// 		}
		// 		pk.add(u);
		// 	}
		// }

		boolean update_volume = true;
		applyWaterEffects(update_volume);
	}
	catch(Exception e) {// singular matrix
		System.out.println("Matrix is singular: reverting to simple brush.");
		pins.clear();
		// brush.deform(P);
	}
}

// void deform(){
// 	// Simply move all points within initial radius of circle over by force
// 	for (PVector p : my_splash.splash){
// 		float sqr_dist = (p.x - brush.position.x)*(p.x - brush.position.x) + (p.y - brush.position.y)*(p.y - brush.position.y);
// 		if(sqr_dist < BRUSH_RADIUS*BRUSH_RADIUS){ 
// 			p.x += brush.force.x;
// 			p.y += brush.force.y;
// 		}
// 	}

// 	if(IMAGE_FLAG){
// 		// Simply move all points within initial radius of circle over by force
// 		for (PVector p : my_splash.splash){
// 			float sqr_dist = (p.x - brish.position.x)*(p.x - brish.position.x) + (p.y - brish.position.y)*(p.y - brish.position.y);
// 			if(sqr_dist < BRUSH_RADIUS*BRUSH_RADIUS){ 
// 				p.x += brish.force.x;
// 				p.y += brish.force.y;
// 			}
// 		}
// 	}

// 	// Should this be it's own method? QUOKKA!
// 	my_splash.mcf(mouseX, mouseY, BRUSH_RADIUS, GEODESIC_FLAG);

// 	// Move onto constraint
// 	my_splash.projectVolumePositions(initial_area);

// 	// Make mesh spacing uniform
// 	my_splash.reParameterize();
// }

// void drawFuture(){
// 	for (PVector p : my_splash.future_splash){
// 		float sqr_dist = (p.x - brush.position.x)*(p.x - brush.position.x) + (p.y - brush.position.y)*(p.y - brush.position.y);
// 		if(sqr_dist < BRUSH_RADIUS*BRUSH_RADIUS){ 
// 			p.x += brush.force.x;
// 			p.y += brush.force.y;
// 		}
// 	}
// }

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




























