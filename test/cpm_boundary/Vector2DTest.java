package cpm_boundary;

import static org.junit.Assert.*;

import org.junit.Test;

public class Vector2DTest {

	@Test
	public void testRotate90_1(){
		Vector2D v = new Vector2D(1,0);
		v.rotate90();
		assertEquals("Returned the wrong x value after rotation",
				0, v.getX());
		assertEquals("Returned the wrong y value after rotation",
				1, v.getY());
	}
	
	@Test
	public void testRotate90_2(){
		Vector2D v = new Vector2D(3,2);
		v.rotate90();
		assertEquals("Returned the wrong x value after rotation",
				-2, v.getX());
		assertEquals("Returned the wrong y value after rotation",
				3, v.getY());
	}
	
	@Test
	public void testRotate270_1(){
		Vector2D v = new Vector2D(1,0);
		v.rotate270();
		assertEquals("Returned the wrong x value after rotation",
				0, v.getX());
		assertEquals("Returned the wrong y value after rotation",
				-1, v.getY());
	}
	
	@Test
	public void testRotate270_2(){
		Vector2D v = new Vector2D(3,2);
		v.rotate270();
		assertEquals("Returned the wrong x value after rotation",
				2, v.getX());
		assertEquals("Returned the wrong y value after rotation",
				-3, v.getY());
	}
	
	@Test
	public void testAdd1(){
		Vector2D v1 = new Vector2D(43,12);
		Vector2D v2 = new Vector2D(90,-234);
		v1.add(v2);
		assertEquals("Returned the wrong x value after addition",
				133, v1.getX());
		assertEquals("Returned the wrong y value after addition",
				-222, v1.getY());
	}

}
