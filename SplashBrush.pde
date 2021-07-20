

public class SplashBrush{

	PVector position;
  float radius;
  float C1, C2, C3;
  float nu = 0.3;
	PVector displacement = new PVector();
  PVector mom_conv_disp = new PVector();
  float mass;
  PVector momentum = new PVector();


  public SplashBrush(float px, float py, float r) 
  {
  	position	= new PVector(px, py, 0.0);
    radius = r;


    C3 = r/(5-6*nu);
    C1 = C3*(3-4*nu);
    C2 = C3*2*(1-nu)*r*r;

    mass = 0;
  }

  void setDisplacementBasedOnNewPosition(PVector new_position){
  	displacement.set(new_position);
    displacement.sub(position);
  }

  // Copy of current u0 value.
  PVector getDisplacementCopy() {
    PVector copy = new PVector();
    copy.set(displacement);
    return copy;
  }

  public void getMomentum(){
    momentum = PVector.mult(displacement, mass);

    println("momentum.x: " + momentum.x + " momentum.y: " + momentum.y);
  }

  void setPosition(PVector new_position){
  	position.set(new_position);
  }

  /** 
   * Displacement at q using the specified uw control
   * 
   * q = point on circle
   * uwControl = force vector
   * 
   */
  public PVector at_withControl(PVector q, PVector uwControl) {
    PVector v    = new PVector(); //sloth "new"
    v.set(q);
    v.sub(position); /// relative position, "x" for now, equally the vector between point on circle and position of brush

    float   radius2   = radius*radius;
    float   radius3   = radius*radius*radius;
    float   r2    = v.magSq();
    float   re2   = r2 + radius2;
    float   invRe = (float)(1./Math.sqrt(re2));
    float   invRe3 = invRe*invRe*invRe;
    float   invRe5 = invRe3*invRe*invRe;

    float diag = C1*invRe + C2*invRe3;
    float C4   = C3*invRe3; 

    // FORCE BLOCK: UU
    PVector u = new PVector();
    acc(u, diag, uwControl);
    acc(u, C4*v.dot(uwControl), v);
    u.z = 0;
    
    // TORQUE BLOCK: WU
    float Cwu = -0.5*((C1+C3)*invRe3 + 3*C2*invRe5);
    u.z = Cwu * (-v.y*uwControl.x + v.x*uwControl.y);
    
    // FORCE BLOCK: UW
    float tau = uwControl.z;
    u.x -= radius3*invRe3* v.y *tau;
    u.y += radius3*invRe3* v.x *tau;
    
    // TORQUE BLOCK: WW
    u.z += 0.5*radius3*invRe5*(3*radius2-re2)*tau;

    //println("u.z: %f", u.z);
    return u;
  }

  /**
   * A(q-p) 3x3 block for 2d displacement and 1d twist.
   * @param q Position to evaluate at.
   * @param A  9-vector to contain A block in row-major order.
   */
  public void getA2D_3x3Block(PVector q, float[] A) 
  {
    if (A.length != 9) throw new IllegalArgumentException("A must be 9-vector");

    // Auu matrix is: (C1/re + C2/re^3)*I + C3/re^3 x x^T

    PVector v    = new PVector(); //sloth "new"
    v.set(q);
    v.sub(position); /// relative position, "x" 

    float   radius2   = radius*radius;
    float   radius3   = radius*radius*radius;
    float   r2    = v.magSq();
    float   re2   = r2 + radius2;
    float   invRe = (float)(1./Math.sqrt(re2));
    float   invRe3 = invRe*invRe*invRe;
    float   invRe5 = invRe3*invRe*invRe;

    // FORCE BLOCK: Auu
    float diag = C1*invRe + C2*invRe3;
    float C4   = C3*invRe3; 
    A[0] = diag + C4*v.x*v.x; // A00
    A[1] =      + C4*v.x*v.y; // A01
    A[3] = A[1];              // A10
    A[4] = diag + C4*v.y*v.y; // A11

    // FORCE BLOCK: Awu
    float Cwu = -0.5*((C1+C3)*invRe3 + 3*C2*invRe5);
    A[6] = -Cwu*v.y;
    A[7] = +Cwu*v.x;

    // TORQUE BLOCK: Auw
    A[2] = radius3*invRe3*(-v.y);
    A[5] = radius3*invRe3*v.x;

    // TORQUE BLOCK: Aww
    A[8] = 0.5*radius3*invRe5*(3*radius2-re2);
  }

  private void acc(PVector sum, float a, PVector v) {
    sum.x += a*v.x;   
    sum.y += a*v.y;   
    sum.z += a*v.z;
  }

}


