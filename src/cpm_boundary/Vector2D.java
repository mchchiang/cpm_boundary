package cpm_boundary;

public class Vector2D {
	
	private int x;
	private int y;
	
	public Vector2D(int x, int y){
		this.x = x;
		this.y = y;
	}
	
	public Vector2D(Vector2D v){
		this.x = v.x;
		this.y = v.y;
	}
	
	public void setX(int x){
		this.x = x;
	}
	public void setY(int y){
		this.y = y;
	}
	
	public int getX(){return x;}
	public int getY(){return y;}
	
	//return the ith component
	public int get(int i){
		if (i == 0){
			return x;
		} else if (i == 1){
			return y;
		} else {
			return 0;
		}
	}

	public void add(Vector2D v){
		this.x += v.x;
		this.y += v.y;
	}
	
	/**
	 * Rotate the vector <i>counterclockwise</i> by 90 degrees for RH systems;
	 * rotate the vector <i>clockwise</i> by 90 degrees for LH systems;
	 */
	public void rotate90(){
		int newX = -this.y;
		int newY = this.x;
		this.x = newX;
		this.y = newY;
	}
	
	/**
	 * Rotate the vector <i>counterclockwise</i> by 270 degrees for RH systems;
	 * rotate the vector <i>clockwise</i> by 270 degrees for LH systems;
	 */
	public void rotate270(){
		int newX = this.y;
		int newY = -this.x;
		this.x = newX;
		this.y = newY;
	}
	
	public boolean equals(Object obj){
		if (obj instanceof Vector2D){
			Vector2D v = (Vector2D) obj;
			return v.x == this.x && v.y == this.y;
		}
		return false;
	}
	
	public String toString(){
		return "( " + this.x + ", " + this.y + " )";
	}
	
	//static methods	
	public static Vector2D add(Vector2D v1, Vector2D v2){
		Vector2D u = new Vector2D(0,0);
		u.x = v1.x + v2.x;
		u.y = v1.y + v2.y;
		return u;
	}
}
