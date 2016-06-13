package cpm_boundary;

public class Vector2D {
	
	private int x;
	private int y;
	
	public Vector2D(int x, int y){
		this.x = x;
		this.y = y;
	}
	
	public int getX(){return x;}
	public int getY(){return y;}
	
	public boolean equals(Object obj){
		if (obj instanceof Vector2D){
			Vector2D v = (Vector2D) obj;
			return v.x == this.x && v.y == this.y;
		}
		return false;
	}
}
