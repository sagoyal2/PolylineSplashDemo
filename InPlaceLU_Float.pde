
/** In Place LU Decomposition for floats.
 * @author Doug James, + Jama proly
 */
public class InPlaceLU_Float {

  /* ------------------------
   Class variables
   * ------------------------ */

  /** Array for internal storage of decomposition.
   @serial internal array storage.
   */
  private float[][] LU;

  /** Row and column dimensions, and pivot sign.
   @serial column dimension.
   @serial row dimension.
   @serial pivot sign.
   */
  private int m, n, pivsign; 

  /** Internal storage of pivot vector.
   @serial pivot vector.
   */
  private int[] piv;

  private boolean haveTransposedLU = false;

  /* ------------------------ Constructor * ------------------------ */

  /** LU Decomposition
   * Uses a "left-looking", dot-product, Crout/Doolittle algorithm.
   *
   *  @param  A   Rectangular matrix (TO BE REFERENCED AND LUD inplace).
   */
  public InPlaceLU_Float(float[][] A) {

    LU  = A;         
    m   = A.length;  

    // Beware 0-by-0 matrices:
    if (m==0)
      n = 0;
    else
      n = A[0].length;

    piv = new int[m];
    for (int i = 0; i < m; i++) {
      piv[i] = i;
    }
    pivsign = 1;
    float[] LUrowi;
    float[] LUcolj = new float[m];

    // Outer loop.
    //System.out.println("DEBUG: START: sum(LU)="+FloatMatrix.sum(LU));

    for (int j = 0; j < n; j++) {

      // Make a copy of the j-th column to localize references.
      for (int i = 0; i < m; i++)   LUcolj[i] = LU[i][j];

      // Apply previous transformations.

      for (int i = 0; i < m; i++) {
        LUrowi = LU[i];

        // Most of the time is spent in the following dot product.

        int kmax = Math.min(i, j);
        float s = 0.f;
        for (int k = 0; k < kmax; k++) {
          s += LUrowi[k]*LUcolj[k];
        }

        LUrowi[j] = LUcolj[i] -= s;
      }

      // Find pivot and exchange if necessary.
      int p = j;
      for (int i = j+1; i < m; i++) {
        float vji = Math.abs(LUcolj[i]);
        float vjp = Math.abs(LUcolj[p]);
        //int pOld = p;
        if (vji > vjp) {
          p = i;
        }
      }
      if (p != j) {
        for (int k = 0; k < n; k++) {
          float t = LU[p][k]; 
          LU[p][k] = LU[j][k]; 
          LU[j][k] = t;
        }
        int k = piv[p]; 
        piv[p] = piv[j]; 
        piv[j] = k;
        pivsign = -pivsign;
      }

      // Compute multipliers.
      if (j < m & LU[j][j] != 0.f) {
        for (int i = j+1; i < m; i++) {
          LU[i][j] /= LU[j][j];
        }
      }
    }

    /// TRANSPOSE LU:
    //if(!haveTransposedLU)  transposeLU();
  }

  /** Transposes if square */
  private void transposeLU() 
  {
    if (m==n) {//Square matrix:
      for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
          float luIJ = LU[i][j];
          LU[i][j] = LU[j][i];
          LU[j][i] = luIJ;
        }
      }
      haveTransposedLU = true;
    }
  }

  /** Square LU matrices can be transposed for faster backsubstition */
  public boolean haveTransposedLU() {
    return haveTransposedLU;
  }

  boolean isNonsingular       = false;
  boolean isNonsingularCalled = false;

  /** Is the matrix nonsingular?
   *  @return     true if U, and hence A, is nonsingular.
   */
  public boolean isNonsingular () {

    for (int j = 0; j < n; j++) {
      if (LU[j][j] == 0.f)
        return false;
    }
    return true;
  }

  /**
   * Only tests if nonSingular once: result of this method is saved 
   * and used in future calls.
   */
  protected boolean isNonsingularQuick() {
    if (isNonsingularCalled) 
      return isNonsingular;
    else {
      for (int j = 0; j < n; j++) {
        if (LU[j][j] == 0.f) {
          isNonsingular = false;
          return false;
        }
      }
      isNonsingular = true;
      isNonsingularCalled = true;///last for thread-safety
      return true;
    }
  }

  //      /** Return lower triangular factor
  //    @return     L
  //      */
  //      public float[][] getL () {

  //    float[][] L = new float[m][n];
  //    for (int i = 0; i < m; i++) {
  //        for (int j = 0; j < n; j++) {
  //      if (i > j) {
  //          L[i][j] = LU[i][j];
  //      } else if (i == j) {
  //          L[i][j] = 1.f;
  //      } else {
  //          L[i][j] = 0.f;
  //      }
  //        }
  //    }
  //    return L;
  //      }

  //      /** Return upper triangular factor
  //    @return     U
  //      */
  //      public float[][] getU () {

  //    float[][] U = new float[n][n];
  //    for (int i = 0; i < n; i++) {
  //        for (int j = 0; j < n; j++) {
  //      if (i <= j) {
  //          U[i][j] = LU[i][j];
  //      } else {
  //          U[i][j] = 0.f;
  //      }
  //        }
  //    }
  //    return U;
  //      }

  /**
   * @return The fragile reference to the underlying LU matrix.
   */
  public float[][] getLUMatrix() {  
    return LU;
  }

  /** Return pivot permutation vector
   @return  A copy of the piv vector.
   */
  public int[] getPivot () {
    int[] p = new int[m];
    for (int i = 0; i < m; i++) {
      p[i] = piv[i];
    }
    return p;
  }


  /** Determinant
   *  @return     det(A)
   *  @exception  IllegalArgumentException  Matrix must be square
   */
  public float det () {
    if (m != n) {
      throw new IllegalArgumentException("Matrix must be square.");
    }
    float d = (float) pivsign;
    for (int j = 0; j < n; j++) {
      d *= LU[j][j];
    }
    return d;
  }

  /** Diagonal of LU matrix */
  public float[] diagonal() {
    if (m != n)
      throw new IllegalArgumentException("Matrix must be square.");

    float[] diag = new float[n];
    for (int j=0; j<n; j++)  diag[j] = LU[j][j];
    return diag;
  }

  /** Solves A*x = b without allocating any memory.
   *  Uses isNonsingularQuick to test matrix singularity.
   *  @param  b   A vector.
   *  @return     x so that L*U*x = b(piv)
   *  @exception  IllegalArgumentException Matrix row dimensions must agree.
   *  @exception  RuntimeException  Matrix is singular.
   */
  public void solve(float[] b, float[] x) {

    if (b.length != m)
      throw new IllegalArgumentException("b has invalid dimension.");
    if (x.length != n)
      throw new IllegalArgumentException("x has invalid dimension.");
    if (!this.isNonsingularQuick())
      throw new RuntimeException("Matrix is singular.");

    // Copy right hand side with pivoting
    for (int i = 0; i < piv.length; i++) 
      x[i] = b[piv[i]]; 

    if (!haveTransposedLU) {/// Standard memory-paging backsub:

      // Solve L*Y = B(piv,:)
      for (int k = 0; k < n; k++) {
        float xk = x[k];
        for (int i = k+1; i < n; i++) {///!!!!paging!!!!
          x[i] -= xk*LU[i][k];
        }
      }
      // Solve U*X = Y;
      for (int k = n-1; k >= 0; k--) {
        x[k] /= LU[k][k];
        float xk = x[k];
        for (int i = 0; i < k; i++) {///!!!!paging!!!!
          x[i] -= xk*LU[i][k];
        }
      }
    } else {/// Transposed fast backsub:

      // Solve L*Y = B(piv,:)
      for (int k = 0; k < n; k++) {
        float xk = x[k];
        for (int i = k+1; i < n; i++) {
          x[i] -= xk*LU[k][i];
        }
      }
      // Solve U*X = Y;
      for (int k = n-1; k >= 0; k--) {
        x[k] /= LU[k][k];
        float xk = x[k];
        for (int i = 0; i < k; i++) {
          x[i] -= xk*LU[k][i];
        }
      }
    }
  }

  /**
   * Solves AX=B and overwrites B with X. 
   * @param B A rectangular matrix consisting of columns to be inverted.
   */
  public void solve(float[][] B) {

    if (B.length != m)
      throw new IllegalArgumentException("b has invalid dimension.");
    if (!this.isNonsingularQuick())
      throw new RuntimeException("Matrix is singular.");

    // Temporary x array:
    float[] x = new float[n];

    for (int column=0; column<B[0].length; column++) {

      // Copy right hand side with pivoting
      for (int i = 0; i < piv.length; i++) 
        x[i] = B[piv[i]][column]; //b[piv[i]]; 

      if (!haveTransposedLU) {/// Standard memory-paging backsub:

        // "Solve L*Y = b{piv}"
        for (int k = 0; k < n; k++) {
          float xk = x[k];
          for (int i = k+1; i < n; i++) {
            x[i] -= xk*LU[i][k];
          }
        }
        // Solve U*X = Y;
        for (int k = n-1; k >= 0; k--) {
          x[k] /= LU[k][k];
          float xk = x[k];
          for (int i = 0; i < k; i++) {
            x[i] -= xk*LU[i][k];
          }
        }
      } else {/// Transposed fast backsub:

        // "Solve L*Y = b{piv}"
        for (int k = 0; k < n; k++) {
          float xk = x[k];
          for (int i = k+1; i < n; i++) {
            x[i] -= xk*LU[k][i];
          }
        }
        // Solve U*X = Y;
        for (int k = n-1; k >= 0; k--) {
          x[k] /= LU[k][k];
          float xk = x[k];
          for (int i = 0; i < k; i++) {
            x[i] -= xk*LU[k][i];
          }
        }
      }

      // Overwrite B[:][column] with x:
      for (int i=0; i<piv.length; i++)   
        B[i][column] = x[i];
    }// column loop
  }

  /**
   * Memory sensitive inverse computation.
   * Only the inverse matrix and two vectors are allocated.
   * The memory efficiency comes at the cost of copying into the vectors.
   */
  public float[][] inverse() {

    float[][] Z  = new float[n][n]; // the inverse
    if (n==0) return Z;
    inverse(Z);
    return Z;
  }

  /**
   * Very memory sensitive inverse generation of the LU decomp'd matrix.
   * Since a matrix is provided by the user, only two vectors are allocated.
   * @param inverse A preallocated array of suitable dimensions to have
   *   the inverse written into.
   */
  public void inverse(float[][] inverse) {

    if (m!=n) throw new Error("m != n --> can't invert matrix");
    if (!isNonsingular()) 
      throw new Error("Matrix was singular --> can't invert.");
    if (inverse.length != n)
      throw new IllegalArgumentException
        ("Provided inverse array has incorrect row length; "+
        "was "+inverse.length+", but should be "+n);
    if (inverse[0].length != n)
      throw new IllegalArgumentException
        ("Provided inverse array has incorrect column length; "+
        "was "+inverse[0].length+", but should be "+n);

    float[][] Z  = inverse;
    float[]   zj = new float[n];
    float[]   ej = new float[n];
    for (int j=0; j<n; j++) {
      // make identity column:
      if (j>0) ej[j-1] = 0.f;
      ej[j] = 1.f;
      // solve:
      solve(ej, zj);
      // copy zj into j'th column of Z:
      for (int i=0; i<n; i++)  Z[i][j] = zj[i];
    }
  }


}