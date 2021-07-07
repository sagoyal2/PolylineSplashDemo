
import java.util.*;

float MESH_THRESHOLD = 10;
float INITIAL_WEIGHT = 0.0002;

public class PolylineSplash{

	ArrayList<PVector> splash;
	private ArrayList<PVector> normals; //unit normals
	private ArrayList<Float> curvatures;
	private ArrayList<PVector> jacobian; //1 by # of mesh points
	ArrayList<Float> weight; //represents alpha_i/mass

	public PolylineSplash(float width, float height, float mesh_resolution, float initial_radius){

		splash 	= new ArrayList<PVector>();
		normals = new ArrayList<PVector>();
		curvatures = new ArrayList<Float>();
		jacobian = new ArrayList<PVector>();
		weight = new ArrayList<Float>();
		for(int i = 0; i < mesh_resolution; i++){
			splash.add(new PVector(	width/2.0+(float)(initial_radius*Math.cos(2*Math.PI*i/(double)mesh_resolution)), 
															height/2.0+(float)(initial_radius*Math.sin(2*Math.PI*i/(double)mesh_resolution)), 
															(float)0));
			weight.add(INITIAL_WEIGHT);
		}

		MESH_THRESHOLD = PVector.dist(splash.get(1), splash.get(2));
	}

	// Copy Constructor
	public PolylineSplash(PolylineSplash copy_polylinesplash){

		splash 	= new ArrayList<PVector>();
		weight = new ArrayList<Float>();

		for(int i = 0; i < copy_polylinesplash.splash.size(); i++){
			splash.add(copy_polylinesplash.splash.get(i).copy());
			weight.add(copy_polylinesplash.weight.get(i)); //I think this is correct?
		}

		normals = new ArrayList<PVector>();
		curvatures = new ArrayList<Float>();
		jacobian = new ArrayList<PVector>();
	}

	public void viewSurface(){
	  fill(0);
	  stroke(0);

	  for (int i = 0; i < splash.size(); i++) {
	    PVector a = splash.get(i);
	    PVector b = splash.get((i+1)%splash.size());
	    line(a.x, a.y, b.x, b.y);
	  }
	}

	public void viewPoints(){
		stroke(0);
		for(int i = 0; i < splash.size(); i++){
			PVector a = splash.get(i);
  		fill(color(52, 125, 235));
  		circle(a.x, a.y, 6);
  	}
	}

	// https://stackoverflow.com/questions/1406029/how-to-calculate-the-volume-of-a-3d-mesh-object-the-surface-of-which-is-made-up/1568551#1568551
	public float getArea(){
		float area = 0.0;
		
		for(int i = 0; i < splash.size(); i++){
			PVector a = splash.get(i);
			PVector b = splash.get((i+1)%splash.size());
			area += -b.x*a.y + a.x*b.y;
  	}
  	
  	return area/2.0;
	}


	// http://chenlab.ece.cornell.edu/Publication/Cha/icip01_Cha.pdf
	private void getNormals(){
		normals.clear();

		// first we calculate surface normal then combine and average neighbors to calculate point normal

		PVector x_axis = new PVector(1.0, 0.0, 0.0);
		PVector y_axis = new PVector(0.0, 1.0, 0.0);

		ArrayList<PVector> surface_normals = new ArrayList<PVector>();
		for(int i = 0; i < splash.size(); i++){
			PVector a = splash.get(i);
			PVector b = splash.get((i+1)%splash.size());

			PVector numerator = PVector.add(PVector.mult(x_axis, -1*(b.y - a.y)),
									 										PVector.mult(y_axis, b.x - a.x));
			float denominator = sqrt(sq(b.x - a.x) + sq(b.y - a.y));
			surface_normals.add(PVector.mult(numerator, denominator));
		}


		for(int i = 0; i < surface_normals.size(); i++){
			PVector a = surface_normals.get((i + surface_normals.size() - 1)%surface_normals.size());
			PVector b = surface_normals.get(i);

			// here we add a -1 because we want the normal to point outwards
			PVector normal_dir = PVector.mult(PVector.add(a, b), -0.5);

			normals.add(normal_dir.normalize());
		}
	}

	public void drawPointNormals(){
		getNormals();

		color normal_color = color(52, 125, 235);
		float normal_scale = 50;

		for(int i = 0; i < normals.size(); i++){
			drawVector(splash.get(i), normals.get(i), normal_color, normal_scale);
			// drawVector(PVector.mult(PVector.add(splash.get(i), splash.get((i+1)%splash.size())), 0.5), normals.get(i), normal_color, normal_scale);
		}
	}


	private void drawVector(PVector position, PVector vector, color c, float scale){
		stroke(c);
		line(position.x, position.y, position.x +  scale*vector.x, position.y + scale*vector.y);
		stroke(0, 0, 0);
	}


	// https://www.youtube.com/watch?v=8JCR6z3GLVI&list=PL9_jI1bdZmz0hIrNCMQW1YmZysAiIYSSS&index=2&ab_channel=KeenanCrane
	// We will maintain the following property: Total curvature flow remains constant throughout the flow
	// other possible: 	Drift (the center of mass does not drift from origin)
	//									Round (flow is stationary for circular curves)

	private void getCurvature(){
		curvatures.clear();

		// Turning Angle:
		// PVector x_axis = new PVector(1.0, 0.0, 0.0);

		// for(int i = 0; i < splash.size(); i++){
		// 	PVector prior 	= splash.get((i + splash.size() - 1)%splash.size());
		// 	PVector current = splash.get(i);
		// 	PVector next 		= splash.get((i + 1)%splash.size());


		// 	PVector leg_one = PVector.sub(current, prior).normalize();
		// 	PVector leg_two = PVector.sub(next, current).normalize();

		// 	float angle_one = acos(PVector.dot(leg_one, x_axis));
		// 	float angle_two = acos(PVector.dot(leg_two, x_axis));

		// 	//here we add a -1 because we want the normal to point outwards, and positive curvature means outwards
		// 	float k = angle_one - angle_two;
		// 	//k*=-1;

		// 	curvatures.add(k);
		// }

		//println("curvature: " + curvatures.get(1));


		// Osculating Circle:
		PVector x_axis = new PVector(1.0, 0.0, 0.0);

		for(int i = 0; i < splash.size(); i++){
			PVector prior 	= splash.get((i + splash.size() - 1)%splash.size());
			PVector current = splash.get(i);
			PVector next 		= splash.get((i + 1)%splash.size());

			PVector leg_one = PVector.sub(current, prior).normalize();
			PVector leg_two = PVector.sub(next, current).normalize();

			float angle_one = acos(PVector.dot(leg_one, x_axis));
			float angle_two = acos(PVector.dot(leg_two, x_axis));

			if(current.y > prior.y) { angle_one *= -1;}
			if(next.y > current.y) { angle_two *= -1;}

			float exterior_angle = angle_two - angle_one;

			float w = PVector.sub(prior, next).mag();

			float k = 2*sin(exterior_angle)/w;
			curvatures.add(k);
		}

	}


	private void drawVectorWithLabel(PVector position, PVector vector, float label, color c, float scale){
		stroke(c);
		fill(c);

		float sign = label/abs(label);

		line(position.x, position.y, position.x +  sign*scale*vector.x, position.y + sign*scale*vector.y);
		text(nf(label,1, 7), position.x + sign*(scale+1)*vector.x, position.y + sign*(scale+1)*vector.y);
		stroke(0, 0, 0);
	}

	public void labelCurvature(){
		getCurvature();
		getNormals();

		color curvature_color = color(237, 164, 245);
		float curvature_scale = 50;

		for(int i = 0; i < splash.size(); i++){
			drawVectorWithLabel(splash.get(i), normals.get(i), curvatures.get(i), curvature_color, curvature_scale);
		}
	}


	public void refineMesh(){
	  int k = 0;

	  while(k < splash.size()){

	  	PVector a = splash.get(k);
	  	PVector b = splash.get((k+1)%splash.size());
	  	float dist = a.dist(b);

	  	//subdivide
	  	if(dist > 2*MESH_THRESHOLD){
	  		PVector mid = new PVector();
	  		mid = PVector.lerp(a, b, 0.5);
	  		splash.add((k+1)%splash.size(), mid);

	  		float mid_weight = 0.5*(weight.get(k) + weight.get((k+1) %weight.size()));
	  		weight.add((k+1)%splash.size(), mid_weight);
	  	}

	  	//coarsen
	  	else if(dist < MESH_THRESHOLD/2.4){
	  		a.lerp(b, 0.5);
	  		splash.remove((k+1)%splash.size());
	  		weight.remove((k+1)%splash.size());
	  		k--;
	  	}

	  	k++;
	  }
	}


	public void showWeight(){
		getNormals();

		color weight_color = color(237, 159, 50);
		float weight_scale = 50;

		for(int i = 0; i < splash.size(); i++){
			drawVectorWithLabel(splash.get(i), normals.get(i), weight.get(i), weight_color, weight_scale);
		}
	}

	public void getJacobian(){
		jacobian.clear();

		jacobian = new ArrayList<PVector>();

		for(int i = 0; i < splash.size(); i++){
			PVector prior 	= splash.get((i + splash.size() - 1)%splash.size());
			// PVector current = splash.get(i);
			PVector next 		= splash.get((i + 1)%splash.size());

			float dC_dx = 0.5*(-1*prior.y + next.y);
			float dC_dy = 0.5*(prior.x + -1*next.x);

			jacobian.add(new PVector(dC_dx, dC_dy, 0.0));
		}
	}

	// Project onto constraint
	public void projectPositions(float initial_area){
		getJacobian();

		// // calculate lambda
		// float current_area = getArea();
		// //println("current_area: " + current_area);


		// float numerator = -1*(current_area - initial_area);

		// float denominator_x = 0.0;
		// float denominator_y = 0.0;
		// for(int i = 0; i < splash.size(); i++){
		// 	denominator_x += jacobian.get(i).x*weight.get(i)*jacobian.get(i).x;
		// 	denominator_y += jacobian.get(i).y*weight.get(i)*jacobian.get(i).y;
		// }

		// float lambda_x = numerator/denominator_x;
		// float lambda_y = numerator/denominator_y;

		// // del position
		// for(int i = 0; i < splash.size(); i++){

		// 	PVector del_position = new PVector(	jacobian.get(i).x*lambda_x*weight.get(i),
		// 																			jacobian.get(i).y*lambda_y*weight.get(i));
		// 	// update position
		// 	splash.get(i).add(PVector.mult(del_position, 1));
		// }


		// calculate lambda
		float current_area = getArea();
		//println("current_area: " + current_area);


		float numerator = -1*(current_area - initial_area);

		float denominator = 0.0;
		for(int i = 0; i < splash.size(); i++){
			denominator += jacobian.get(i).x*weight.get(i)*jacobian.get(i).x;
			denominator += jacobian.get(i).y*weight.get(i)*jacobian.get(i).y;
		}

		float lambda = numerator/denominator;

		// del position
		for(int i = 0; i < splash.size(); i++){

			PVector del_position = new PVector(	jacobian.get(i).x*lambda*weight.get(i),
																					jacobian.get(i).y*lambda*weight.get(i));
			// update position
			splash.get(i).add(PVector.mult(del_position, 0.001));
		}

	}
} 































































