
/**
 * Amortized solver for a pin-constrained brush. Supports position and orientation constraints. 
 * Low-rank updates are used to support fast solves.
 * 
 * @author Doug James, Pixar, April 2016.
 */
class FastPinConstraintSolver2D 
{
  private ArrayList<SplashBrush> brushes = new ArrayList<SplashBrush>();
  private float[] c;// controls (fx, fy, tz, ...)

  /** 
   * Initializes solver to a single unconstrained brush.
   */
  public FastPinConstraintSolver2D(SplashBrush brush) 
  {
    brushes.add(brush);
  }

  public void addPin(SplashBrush pin) {
    brushes.add(pin);
    c = null;
  }


  public FastPinConstraintSolver2D(){}

  public void addRig(SplashBrush rig){
    brushes.add(rig);
  }

  public void solve(boolean with_rigs) 
  {
    //long t0 = -System.nanoTime();

    /// SLOTH: VANILLA SOLVER FOR NOW:

    // COLLECT BRUSHES:
    SplashBrush[] B = (SplashBrush[])brushes.toArray(new SplashBrush[0]);

    /// SETUP A c = u = (ux0 uy0 wz0 0 ... 0)':
    int n = B.length;
    int N = 3*n;
    float[][] A = new float[N][N];
    c = new float[N];
    float[]   u = new float[N];

    // COLLECT POSITIONS:
    PVector[] p = new PVector[n];
    for (int k=0; k<n; k++) 
      p[k] = B[k].position;

    // MATRIX:
    float[] aij = new float[9];
    for (int i=0; i<n; i++) { 
      for (int j=0; j<n; j++) {
        B[j].getA2D_3x3Block(p[i], aij);
        for (int i3=0; i3<3; i3++) {
          for (int j3=0; j3<3; j3++) {
            A[3*i+i3][3*j+j3] = aij[3*i3+j3];
          }
        }
      }
    }
    //t0 += System.nanoTime();
    //print("Assembly: "+t0/1000.+"ms;     ");

    //t0 = -System.nanoTime();
    // RHS: 
    
    if(with_rigs){
      for(int i = 0; i < n; i++){
        SplashBrush curr_rig = brushes.get(i);
        PVector u0 = curr_rig.getDisplacementCopy();
        u[i*3 + 0] = u0.x;
        u[i*3 + 1] = u0.y;
        u[i*3 + 2] = u0.z;
      }
    }
    else{
      PVector u0 = brush.getDisplacementCopy();
      u[0] = u0.x;  
      u[1] = u0.y;
      u[2] = u0.z; // wz control (natural twist-free brush --> 0)
    }

    // SOLVE:
    InPlaceLU_Float LU = new InPlaceLU_Float(A);
    LU.solve(u, c);//u=RHS, c=SOLUTION
  }

  /** Deforms point using current control variables. Assumes brushes parameters do not change. */
  public void deform(PVector point) 
  {
    PVector u = getDeformDisplacement(point);
    point.add(u);//Deform point
  }
  /** Returns displacement for specified point using current control variables. Assumes brushes parameters do not change. */
  public PVector getDeformDisplacement(PVector point) 
  {
    if (c==null & brushes.size()>1) throw new RuntimeException("Solve system first.");

    /// Accumulate point's displacements into "u" due to each brush using "c" controls:
    PVector u         = new PVector();
    PVector uwControl = new PVector();
    for (int j=0; j<brushes.size(); j++) {//for each brush...
      uwControl.set(c[3*j], c[3*j+1], c[3*j+2]);

      //u.add(uwControl)

      u.add ( brushes.get(j).at_withControl(point, uwControl) );
    }
    return u;
  }
}




