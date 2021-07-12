









/** 
 *  Based off of:
 * 	Particle-Based Fluid Simulation for Interactive Applications
 *  Link: https://matthias-research.github.io/pages/publications/sca03.pdf
 * 
 *  SplashBrushDemo
 *  @author Samaksh (Avi) Goyal, 6/27/2021
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 	PolylineSplash Class:
 * 		- getArea()
 * 		- refine()
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 */

float BRUSH_RADIUS = 50.0;
float MESH_RESOLUTION = 100;
float INITIAL_DROPLET_RADIUS = 300;

boolean NORMAL_FLAG = false;
boolean CURVATURE_FLAG = false;
boolean WEIGHT_FLAG = false;
boolean WEIGHT_MODE = false;
boolean MCF_FLAG = false;
int mcf_clock = 100;
boolean GEODESIC_FLAG = false;
boolean IMAGE_FLAG = false;


float initial_area = -1.0;
SplashBrush brush;
SplashBrush brish;
PolylineSplash my_splash;
LinkedList<PolylineSplash> undo_splash;

// Drivers
/////////////////////////////////////////////////////////////////
void setup(){
  size(1280, 1024);  
  smooth(8);
  reset();
}

void reset(){ 

	// take in shape parameter?
	my_splash 		= new PolylineSplash(width, height, MESH_RESOLUTION, INITIAL_DROPLET_RADIUS);
	undo_splash 	= new LinkedList<PolylineSplash>();

	initial_area = my_splash.getArea();
  println("initial_area: " + initial_area);
}


void draw(){
  background(255);

  fill(0);
  text("CONTROLS: Radius[w-s], undo[z], show_normals[n], label_curvature[k], weight[w], adjust weight[a], MCF [m], geodesic[g]", 5, 10); 
  text("#UNDOS: " + undo_splash.size(), 5, 25);


  // Visualize
  drawBrushes();

  my_splash.refineMesh();
  my_splash.viewPoints();
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
  	if(mcf_clock < 250){
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

  	// Make mesh spacing uniform
  	my_splash.reParameterize();

  	// Move onto constraint
  	my_splash.projectPositions(initial_area);
  }
  if(GEODESIC_FLAG){
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
  fill(0);

  // float area = my_splash.getArea();
  // println("Area: " + area);

 //  // Time Step with forces
 //  calculate_density();
 //  calculate_color_field(); //not actual color but rather 1/0

 //  zero_force_buffer();

	// set_force_pressure();
 //  set_force_surface_tension();
  // set_force_external();
  // check_boundary();

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
  // if (key == 'w'){
  //   BRUSH_RADIUS *= 1.1;
  //   // brush.setRadius(eps);
  // }
  // if (key == 's'){ 
  //   BRUSH_RADIUS *= 0.9; 
  //   // brush.setRadius(eps);
  // }
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
  	if(IMAGE_FLAG){
  		// brish = new SplashBrush(mouseX, height/2.0 + mouseY);
  	}
  }
}

void mouseWheel(MouseEvent event) {

	float e = event.getCount();
	if (WEIGHT_MODE){
		float weigth_scale;

		if(e > 0.0){
			weigth_scale = 1.1;
		}else{
			weigth_scale = 0.9; 
		}

		reweight(weigth_scale);
	}
	else{
		if(e > 0.0){
			BRUSH_RADIUS *= 1.1;
		}else{
			BRUSH_RADIUS *= 0.9; 
		}
	}
}


/////////////////////////////////////////////////////////////////
void drawBrushes() {
  drawBrush(mouseX, mouseY);

  if(IMAGE_FLAG){
    drawBrush(mouseX, height - mouseY);
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

  if(GEODESIC_FLAG){
  	my_splash.showGeodesic(x, y, BRUSH_RADIUS);
  }
}

void mousePressed() {
  brush = new SplashBrush(mouseX, mouseY);
  
  if(IMAGE_FLAG){
	  brish = new SplashBrush(mouseX, height - mouseY);
  }

  saveSplashState();
}

void saveSplashState() {
  PolylineSplash previous = new PolylineSplash(my_splash);
  undo_splash.addFirst(previous);
}

void mouseDragged() {
  PVector new_position = new PVector(mouseX, mouseY, 0);
  brush.setForceBasedOnNewPosition(new_position);

  PVector new_position_2 =  new PVector(mouseX, height - mouseY, 0);
  if(IMAGE_FLAG){
  	brish.setForceBasedOnNewPosition(new_position_2);
  }

  // solveAndDeform();
  // move points within brush somehow
  deform();

  // Slide brush to newP:
  brush.setPosition(new_position);

  if(IMAGE_FLAG){
  	brish.setPosition(new_position_2);
  }
}

void deform(){
	// Simply move all points within initial radius of circle over by force
	for (PVector p : my_splash.splash){
		float sqr_dist = (p.x - brush.position.x)*(p.x - brush.position.x) + (p.y - brush.position.y)*(p.y - brush.position.y);
		if(sqr_dist < BRUSH_RADIUS*BRUSH_RADIUS){ 
			p.x += brush.force.x;
			p.y += brush.force.y;
		}
  }

  if(IMAGE_FLAG){
  	// Simply move all points within initial radius of circle over by force
		for (PVector p : my_splash.splash){
			float sqr_dist = (p.x - brish.position.x)*(p.x - brish.position.x) + (p.y - brish.position.y)*(p.y - brish.position.y);
			if(sqr_dist < BRUSH_RADIUS*BRUSH_RADIUS){ 
				p.x += brish.force.x;
				p.y += brish.force.y;
			}
	  }
  }

  // Move onto constraint
  my_splash.projectPositions(initial_area);

  // Make mesh spacing uniform
  my_splash.reParameterize();
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



// void draw() {
//   PVector dir = PVector.sub(mx(), pmx());
//   float magnitude = dir.mag();
//   background(255);
//   fill(0);
//   String info = "Degrees: " + int(degrees(dir.heading())) + "\nMagnitude: " + int(magnitude);
//   text(info, 5, 5);
//   translate(width/2, height/2);
//   line(0, 0, dir.x, dir.y);
// }
 
// PVector mx() {
//   return new PVector(mouseX, mouseY);
// }
 
// PVector pmx() {
//   return new PVector(pmouseX, pmouseY);
// }






















