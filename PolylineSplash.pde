
import java.util.*;

float MESH_THRESHOLD = 10;
float INITIAL_WEIGHT = 50;

public class PolylineSplash{

	ArrayList<PVector> splash;
	ArrayList <PVector> future_splash;
	ArrayList<PVector> normals; //unit normals
	private ArrayList<Float> curvatures;
	private ArrayList<PVector> jacobian; //1 by # of mesh points
	ArrayList<Float> weight; //represents alpha_i/mass
	ArrayList<Float> geodesic;
	ArrayList<Float> depth;

	public PolylineSplash(float width, float height, float mesh_resolution, float initial_radius){

		splash 	= new ArrayList<PVector>();
		future_splash 	= new ArrayList<PVector>();
		normals = new ArrayList<PVector>();
		curvatures = new ArrayList<Float>();
		jacobian = new ArrayList<PVector>();
		weight = new ArrayList<Float>();
		geodesic = new ArrayList<Float>();
		depth = new ArrayList<Float>();

		float center_x = 7*width/8.0;
		float center_y = 7*height/8.0;

		for(int i = 0; i < mesh_resolution; i++){
			splash.add(new PVector(	center_x+(float)(initial_radius*Math.cos(2*Math.PI*i/(double)mesh_resolution)), 
															center_y+(float)(initial_radius*Math.sin(2*Math.PI*i/(double)mesh_resolution)), 
															(float)0));
			future_splash.add(splash.get(i).copy());
			weight.add(INITIAL_WEIGHT);
		}

		MESH_THRESHOLD = PVector.dist(splash.get(1), splash.get(2));
	}

	// Copy Constructor
	public PolylineSplash(PolylineSplash copy_polylinesplash){

		splash 	= new ArrayList<PVector>();
		future_splash 	= new ArrayList<PVector>();
		weight = new ArrayList<Float>();

		for(int i = 0; i < copy_polylinesplash.splash.size(); i++){
			splash.add(copy_polylinesplash.splash.get(i).copy());
			future_splash.add(copy_polylinesplash.future_splash.get(i).copy());
			weight.add(copy_polylinesplash.weight.get(i)); //I think this is correct?
		}

		normals = new ArrayList<PVector>();
		curvatures = new ArrayList<Float>();
		jacobian = new ArrayList<PVector>();
		geodesic = new ArrayList<Float>();
		depth = new ArrayList<Float>();
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

	public void viewPoints(boolean SHOW_FUTURE, boolean WITH_MOM_CONSTRAINT){

		color blue = color(0, 145 , 255); //bigger 		-- incompressable
		color red = color(255, 68, 0);		//smaller		-- compressable

		stroke(0);
		for(int i = 0; i < splash.size(); i++){
			PVector a = splash.get(i);
			float scaling = (1-(weight.get(i))/100);
			fill(lerpColor(blue, red, scaling));
			circle(a.x, a.y, 6);

			if(SHOW_FUTURE){
				PVector b = future_splash.get(i);
				
				if(WITH_MOM_CONSTRAINT){
					fill(52, 174, 235);
				}
				else{
					fill(255, 128, 89);
				}
				circle(b.x, b.y, 6);
			
				// fill(232, 233, 237);
				// if(PVector.dist(a,b) > 0.0){
				// 	line(a.x, a.y, b.x, b.y);
				// }
			}
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

	//https://stackoverflow.com/questions/217578/how-can-i-determine-whether-a-2d-point-is-within-a-polygon
	public boolean getIntersections(PVector source){

		int intersection_count = 0;

		for(int j = 0; j < splash.size(); j++){

			PVector a = splash.get(j);
			PVector b = splash.get((j+1)%splash.size());

			float m1 = (b.y - a.y)/(b.x - a.x);
			float m2 = source.y/source.x;

			float x_sol = (m1*a.x - m2*source.x - a.y + source.y)/(m1 - m2);
			float y_sol = m1*(x_sol - a.x) + a.y;

			PVector target = new PVector(x_sol, y_sol, 0.0);
			// float dist = PVector.dist(source, target);

			// PVector result = PVector.sub(target, source).normalize();
			// float angle = acos(PVector.dot(n, result));
			// && (angle > 0)
			if((abs((PVector.dist(a, target) + PVector.dist(target, b)) - PVector.dist(a, b)) < 1) ){
				intersection_count++;
			}
		}

		// Early Out
		if((intersection_count == 0) || (intersection_count%2 ==1)){
			return false;
		}


		//https://stackoverflow.com/questions/10673740/how-to-check-if-a-point-x-y-is-inside-a-polygon-in-the-cartesian-coordinate-sy
		float angle = 0;
		PVector p1, p2;

		for(int i=0 ; i<splash.size();i++){
			PVector a = splash.get(i);
			PVector b = splash.get((i+1)%splash.size());
			p1 = PVector.sub(a, source);
			p2 = PVector.sub(b, source);

			float theta1 = atan2(p1.y, p1.x);
			float theta2 = atan2(p2.y, p2.x);

			float dtheta = theta2 - theta1;
			if(dtheta > PI){
				dtheta -= 2*PI;
			}
			if(dtheta < -1.*PI){
				dtheta += 2*PI;
			}

			angle += dtheta;
		}

		if(abs(angle) < PI){
			return false;
		}else{
			return true;
		}

	}

	public void showVolume(){
		//for each pixel color it blue if inside

		color blue = color(52, 125, 235);
		for(int i = 1; i < width; i++){
			for(int j = 1; j < height; j++){

				boolean hit = getIntersections(new PVector(i, j, 0.0));
				if(hit){
					set(i, j, blue);
				}
			}
		}
	}

	// http://chenlab.ece.cornell.edu/Publication/Cha/icip01_Cha.pdf
	public void getNormals(){
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
			// drawVector(splash.get(i), normals.get(i), normal_color, normal_scale);
			// drawVector(PVector.mult(PVector.add(splash.get(i), splash.get((i+1)%splash.size())), 0.5), normals.get(i), normal_color, normal_scale);
			drawVectorWithLabel(splash.get(i), normals.get(i), i, normal_color, normal_scale);
		}
	}


	private void drawVector(PVector position, PVector vector, color c, float scale){
		stroke(c);
		line(position.x, position.y, position.x +  scale*vector.x, position.y + scale*vector.y);
		stroke(0, 0, 0);
	}

	private void drawVectorWithLabel(PVector position, PVector vector, float label, color c, float scale){
		stroke(c);
		fill(c);

		float sign = label/abs(label);

		line(position.x, position.y, position.x +  sign*scale*vector.x, position.y + sign*scale*vector.y);
		// text(nf(label, 2, 2), position.x + sign*(scale+1)*vector.x, position.y + sign*(scale+1)*vector.y);
		text(label,  position.x + sign*(scale+1)*vector.x, position.y + sign*(scale+1)*vector.y);
		stroke(0, 0, 0);
	}

	// https://www.youtube.com/watch?v=8JCR6z3GLVI&list=PL9_jI1bdZmz0hIrNCMQW1YmZysAiIYSSS&index=2&ab_channel=KeenanCrane
	// We will maintain the following property: Total curvature flow remains constant throughout the flow
	// other possible: 	Drift (the center of mass does not drift from origin)
	//									Round (flow is stationary for circular curves)

	private void getCurvature(){
		curvatures.clear();

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
			if(dist > 1.1*MESH_THRESHOLD){
				PVector mid = new PVector();
				mid = PVector.lerp(a, b, 0.5);
				splash.add((k+1)%splash.size(), mid);
				future_splash.add((k+1)%splash.size(), mid);

				float mid_weight = 0.5*(weight.get(k) + weight.get((k+1) %weight.size()));
				weight.add((k+1)%splash.size(), mid_weight);
			}

			//coarsen
			else if(dist < MESH_THRESHOLD/2.4){
				a.lerp(b, 0.5);
				splash.remove((k+1)%splash.size());
				future_splash.remove((k+1)%splash.size());
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

	private void getVolumeJacobian(){
		jacobian.clear();

		jacobian = new ArrayList<PVector>();

		for(int i = 0; i < splash.size(); i++){
			PVector prior 	= splash.get((i + splash.size() - 1)%splash.size());
			// PVector current = splash.get(i);
			PVector next 		= splash.get((i + 1)%splash.size());

			float dC_dx = 0.5*(-1.0*prior.y + next.y);
			float dC_dy = 0.5*(prior.x + -1.0*next.x);

			jacobian.add(new PVector(dC_dx, dC_dy, 0.0));
		}
	}

	// Project onto constraint
	public void projectVolumePositions(float initial_area){
		getVolumeJacobian();

		// calculate lambda
		float current_area = getArea();
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
			splash.get(i).add(PVector.mult(del_position, 1));
		}
	}

	public void mcf(float x, float y, float brush_radius, boolean g_flag){
		getNormals();
		getCurvature();
		getGeodesic(x, y, brush_radius);

		float flow_scale = 10;
		
		// Explict Euler
		for(int i = 0; i < splash.size(); i++){
			float k = curvatures.get(i);
			PVector normal_dir = normals.get(i);

			float geo_scale = 1;
			if(g_flag){
				geo_scale *= geodesic.get(i);
			}

			splash.get(i).add(PVector.mult(normal_dir, geo_scale*flow_scale*k));
		}
	}

	public void getGeodesic(float x, float y, float brush_radius){
		geodesic.clear();

		PVector center = new PVector(x, y, 0.0);
		int closest_point_index = indexOfClosestPoint(center);

		for (int k=0; k < splash.size(); k++)	geodesic.add(Float.MAX_VALUE);
		geodesic.set(closest_point_index, 0.0);
		
		// SWEEP DOWN/CCW:
		float distSum = 0;
		PVector pLast = splash.get(closest_point_index);
		int     iNext = closest_point_index;
		for (int k=1; k < geodesic.size(); k++) {
			iNext--;
			if (iNext < 0) iNext = geodesic.size()-1;
			PVector pNext = splash.get(iNext);
			distSum += pNext.dist(pLast);
			if (geodesic.get(iNext) > distSum)  geodesic.set(iNext, distSum);
			pLast = pNext;
		}

		// Sweep UP/CW;
		distSum = 0;
		pLast		=	splash.get(closest_point_index);
		iNext 	= closest_point_index;
		for (int k=1; k < geodesic.size(); k++) {
			iNext++;
			if (iNext >= geodesic.size()) iNext = 0;
			PVector pNext = splash.get(iNext);
			distSum += pNext.dist(pLast);
			if (geodesic.get(iNext) > distSum)  geodesic.set(iNext, distSum);
			pLast = pNext;
		}	


		// if within radius then set distance to 0;
		for (int k = 0; k < geodesic.size(); k++){
			if(center.dist(splash.get(k)) < brush_radius){
				geodesic.set(k, 0.0);
			}
		}		

		// scale geodesic weight
		float maximum_dist = Float.MIN_VALUE;
		for (int k = 0; k < geodesic.size(); k++){
			if(geodesic.get(k) > maximum_dist){
				maximum_dist = geodesic.get(k);
			}
		}

		for (int k = 0; k < geodesic.size(); k++){
			float g = geodesic.get(k);
			geodesic.set(k, pow((maximum_dist - g)/(maximum_dist), 5));
		}		
	}

	public int indexOfClosestPoint(PVector center) {
		int   minIndex = -1;
		float minDist  = Float.MAX_VALUE;
		for (int k=0; k<splash.size(); k++) {
			PVector pk = splash.get(k);
			float   dk = center.dist(pk);//sloth sqrt
			if (dk < minDist) {
				minDist = dk;
				minIndex = k;
			}
		}
		return minIndex;
	}

	public void showGeodesic(float x, float y, float brush_radius){
		getGeodesic(x, y, brush_radius);
		getNormals();

		float generic_scale = 50;
		color geo_color = color(46, 166, 50);

		for(int i = 0; i < splash.size(); i++){
			float geo_scale = generic_scale * geodesic.get(i);
			drawVector(splash.get(i), normals.get(i), geo_color, geo_scale);
		}
	}

	// Make mesh spacing uniform
	public void reParameterize(){

		int iteration = 1;
		while(iteration > 0){
			// Move points to be average of neighbors
			int k = 0;
			while(k < splash.size()){
				PVector prior 	= splash.get((k + splash.size() - 1)%splash.size());
				PVector current = splash.get(k);
				PVector next 		= splash.get((k + 1)%splash.size());

				PVector centroid = new PVector();
				centroid = PVector.add(PVector.add(PVector.mult(current, 1.0), prior), next);
				current.set(PVector.mult(centroid, 1.0/3.0));
				k++;
			}
			iteration--;
		}
	}

	public void startFutureSplash(){
		future_splash.clear();

		for(PVector p: splash){
			future_splash.add(p.copy());
		}
	}

	public void projectToFuture(){
		// modified geodesic
		geodesic.clear();
		int closest_point_index = 0;
		for(int i = 0; i < splash.size(); i++){
			PVector a = splash.get(i);
			PVector b = future_splash.get(i);
			if(PVector.dist(a,b) > 0.0){
				closest_point_index = i;
				geodesic.add(0.0);
			}else{
				geodesic.add(Float.MAX_VALUE);
			}
		}
		
		// SWEEP DOWN/CCW:
		float distSum = 0;
		PVector pLast = splash.get(closest_point_index);
		int     iNext = closest_point_index;
		for (int k=1; k < geodesic.size(); k++) {
			if (Math.signum(geodesic.get(iNext)) == 0){ distSum = 0;}
			iNext--;
			if (iNext < 0) iNext = geodesic.size()-1;
			PVector pNext = splash.get(iNext);
			distSum += pNext.dist(pLast);
			if (geodesic.get(iNext) > distSum)  geodesic.set(iNext, distSum);
			pLast = pNext;
		}

		// Sweep UP/CW;
		distSum = 0;
		pLast		=	splash.get(closest_point_index);
		iNext 	= closest_point_index;
		for (int k=1; k < geodesic.size(); k++) {
			if (Math.signum(geodesic.get(iNext)) == 0){ distSum = 0;}
			iNext++;
			if (iNext >= geodesic.size()) iNext = 0;
			PVector pNext = splash.get(iNext);
			distSum += pNext.dist(pLast);
			if (geodesic.get(iNext) > distSum)  geodesic.set(iNext, distSum);
			pLast = pNext;
		}	
	
		// scale geodesic weight
		float maximum_dist = Float.MIN_VALUE;
		for (int k = 0; k < geodesic.size(); k++){
			if(geodesic.get(k) > maximum_dist){
				maximum_dist = geodesic.get(k);
			}
		}

		// set weight based on geodesic
		for (int k = 0; k < geodesic.size(); k++){
			float g = geodesic.get(k);
			weight.set(k, 100*pow((maximum_dist - g)/(maximum_dist), 2));
		}	
		////////////////////////////////////////////////////

		// move all current splash points to their new location
		for(int i = 0; i < splash.size(); i++){
			splash.set(i, future_splash.get(i).copy());
		}
	}

	public void setWeightToNeutral(){
		for(int i = 0; i < weight.size(); i++){
			weight.set(i, INITIAL_WEIGHT);
		}
	}


	void getDepth(){
		depth.clear();
		getNormals();

		// for(int i = 0; i < splash.size(); i++){
		// 	PVector source = splash.get(i);
		// 	PVector n = normals.get(i);

		// 	float min_depth = Float.MAX_VALUE;

		// 	for(int j = 0; j < splash.size(); j++){
		// 		if((i == j) || (i == (j+1)%splash.size())){
		// 			continue;
		// 		}

		// 		PVector a = splash.get(j);
		// 		PVector b = splash.get((j+1)%splash.size());

		// 		float m1 = (b.y - a.y)/(b.x - a.x);
		// 		float m2 = n.y/n.x;

		// 		float x_sol = (m1*a.x - m2*source.x - a.y + source.y)/(m1 - m2);
		// 		float y_sol = m1*(x_sol - a.x) + a.y;

		// 		PVector target = new PVector(x_sol, y_sol, 0.0);
		// 		float dist = PVector.dist(source, target);

		// 		PVector result = PVector.sub(target, source).normalize();
		// 		float angle = acos(PVector.dot(n, result));
				
		// 		// && (angle > 0)
		// 		if((abs(PVector.dist(a, target) + PVector.dist(target, b) - PVector.dist(a, b)) < 5) ){
		// 			min_depth = dist;
		// 		}
		// 	}

		// 	if(min_depth > 1024){
		// 		min_depth = -1;
		// 	}
		// 	depth.add(min_depth);
		// }


		for(int i = 0; i < splash.size(); i++){
			PVector source = splash.get(i);
			PVector n = normals.get(i);


			//float min_angle = PI;
			float min_depth = Float.MAX_VALUE;

			for(int j = 0; j < splash.size(); j++){
				if(i == j){continue;}

				PVector target = splash.get(j);
				PVector result = PVector.sub(target, source).normalize();
				float angle = acos(PVector.dot(PVector.mult(n, -1.0), result));
				float depth = PVector.dist(target, source);

				if((abs(angle) < 0.5) && depth < min_depth){
					//min_angle = angle;

					min_depth = depth;
				}
			}

			depth.add(min_depth);
		}
	}

	public void showDepth(){
		getNormals();
		getDepth();

		color depth_color = color(252, 132, 3);
		float depth_scale = 50;

		for(int i = 0; i < splash.size(); i++){
			// if(i == 0){
			// drawVectorWithLabel(splash.get(i), normals.get(i), depth.get(i), depth_color, -1*depth.get(i));

			// }
			drawVectorWithLabel(splash.get(i), normals.get(i), depth.get(i), depth_color, -1*depth.get(i));
		}
	}

	public void fixDepth(){
		getDepth();

		boolean inversion_present = true;

		// QUOKKA - perhaps don't make this a while loop, and have it do something better
		//while(inversion_present){
			//inversion_present = false;
			for(int i = 0; i < splash.size(); i++){
				PVector p = splash.get(i);
				if(depth.get(i) < MESH_THRESHOLD){
					inversion_present = true;
					p.x += 5*normals.get(i).x;
					p.y += 5*normals.get(i).y;
				}
			}
		//}
	}



	public PVector findNecks(){
		float NECK_THRESHOLD = 80;
		Set<PVector> points_to_cluster = new HashSet<PVector>();

		for(int i = 0; i < depth.size(); i++){
			if(depth.get(i) < NECK_THRESHOLD){
				points_to_cluster.add(splash.get(i));
			}
		}

		//find neck
		PVector neck = new PVector();
		for(PVector p : points_to_cluster){
			neck.add(p);
		}

		neck.div(points_to_cluster.size());
		return neck;

		// right way to do this is to find number of clusters, and then group points in said clusters
	}



} 






















































