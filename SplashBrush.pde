

public class SplashBrush{

	PVector position;
	PVector force = new PVector();


  public SplashBrush(float px, float py) 
  {
  	position	= new PVector(px, py, 0.0);
  }

  void setForceBasedOnNewPosition(PVector new_position){
  	force.set(new_position);
    force.sub(position);
  }

  void setPosition(PVector new_position){
  	position.set(new_position);
  }


}