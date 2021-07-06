









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
float MESH_RESOLUTION = 200;
float INITIAL_DROPLET_RADIUS = 300;
boolean NORMAL_FLAG = false;
boolean CURVATURE_FLAG = false;


SplashBrush brush;
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
}


void draw(){
  background(255);

  fill(0);
  text("CONTROLS: Radius[w-s], undo[z], show_normals[n], label_curvature[k]", 5, 10); 
  text("#UNDOS: " + undo_splash.size(), 5, 25);


  // Visualize
  drawBrushes();

  //PolylineSplash my_splash = new PolylineSplash(width, height, MESH_RESOLUTION, INITIAL_DROPLET_RADIUS);

  my_splash.viewPoints();
  if(NORMAL_FLAG){
  	my_splash.drawPointNormals();
  }
  if(CURVATURE_FLAG){
  	my_splash.labelCurvature();
  }

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
  if (key == 'w'){
    BRUSH_RADIUS *= 1.1;
    // brush.setRadius(eps);
  }
  if (key == 's'){ 
    BRUSH_RADIUS *= 0.9; 
    // brush.setRadius(eps);
  }
  if (key == 'z'){//undo
    if (undo_splash.size() > 0) {
      my_splash = undo_splash.pollFirst();
    }
  }
}


/////////////////////////////////////////////////////////////////
void drawBrushes() {
  drawBrush(mouseX, mouseY);
}

void drawBrush(float x, float y)
{
  noFill();
  circle(x, y, 2*BRUSH_RADIUS);
}

void mousePressed() {
  brush = new SplashBrush(mouseX, mouseY);
  saveSplashState();
}

void saveSplashState() {
  PolylineSplash previous = new PolylineSplash(my_splash);
  undo_splash.addFirst(previous);
}

void mouseDragged() {
  PVector new_position = new PVector(mouseX, mouseY, 0);
  brush.setForceBasedOnNewPosition(new_position);

  // solveAndDeform();
  // move points within brush somehow
  deform();

  // Slide brush to newP:
  brush.setPosition(new_position);
}


void deform(){
	//simply move all points within initial radius of circle over by force
	for (PVector p : my_splash.splash){
		float sqr_dist = (p.x - brush.position.x)*(p.x - brush.position.x) + (p.y - brush.position.y)*(p.y - brush.position.y);
		if(sqr_dist < BRUSH_RADIUS*BRUSH_RADIUS){ 
			p.x += brush.force.x;
			p.y += brush.force.y;
		}
  }
}



















































