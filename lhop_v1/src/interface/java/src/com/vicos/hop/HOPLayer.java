package com.vicos.hop;

import java.util.Iterator;

public class HOPLayer {

	public static class HOPLinkIterator {
		
		private HOPLayer top, bottom;
		
		public HOPLinkIterator(HOPLayer top) {
			this.bottom = top.parent;
			this.top = top;
			
			while (top.links[id] == null) {
				id ++;
				if (id >= top.size())
					return;
			}
			
			linkedId = top.links[id][link];
			
		}
		
		private int id = 0, link = 0, linkedId;
		
		public int getIndex() {
			return id;
		}
		
		public boolean next() {
			
			if (id >= top.size())
				return false;
			
			if (link < top.links[id].length - 1) {
				
				link++;
				linkedId = top.links[id][link];
				
				return true;
				
			} else {
				id++;
				if (id >= top.size())
					return false;
				while (top.links[id] == null) {
					id ++;
					if (id >= top.size())
						return false;
				}
				
				link = 0;
				linkedId = top.links[id][link];
				return true;
			}

			//return true;
		}
		
		public double getTopX() {
			return top.x[id];
		}
		
		public double getTopY() {
			return top.y[id];
		}
		
		public double getTopWeight() {
			return top.w[id];
		}
		
		public double getBottomX() {
			return bottom.x[linkedId];
		}
		
		public double getBottomY() {
			return bottom.y[linkedId];
		}
		
		public double getBottomWeight() {
			return bottom.w[linkedId];
		}
	}
	
	public class HOPPartIterator {
		
		private int id = 0;
		
		public boolean next() {
			id ++;
			if (id >= size())
				return false;

			return true;
		}
		
		public double getX() {
			return x[id];
		}
		
		public double getY() {
			return y[id];
		}
		
		public double getWeight() {
			return w[id];
		}
		
		public int getLibraryReference() {
			return references[id];
		}
		
		public int getIndex() {
			return id;
		}
	}
	
	public class HOPPart implements Iterable<HOPPart> {
		
		private int id;
		
		protected HOPPart(int id) {
			this.id = id;
		}
		
		public double getX() {
			return x[id];
		}
		
		public double getY() {
			return y[id];
		}
		
		public double getWeight() {
			return w[id];
		}
		
		public int getLibraryReference() {
			return references[id];
		}

		@Override
		public Iterator<HOPPart> iterator() {
			return new Iterator<HOPPart>() {
				
				private int current = 0;
				
				@Override
				public void remove() {}
				
				@Override
				public HOPPart next() {
					if (!hasNext())
						return null;
				
					int lid = links[id][current];
					
					return parent.get(lid);
					
				}
				
				@Override
				public boolean hasNext() {
					if (parent == null)
						return false;
					
					return (current < links[id].length);
					
				}
			};
		}
		
	}

	private double[] x, y, w;
	
	private int[] references;
	
	private int[][] links;
	
	private HOPLayer parent = null;

	private void allocateLinks(int part, int size) {
		
		links[part] = new int[size];
		
	}
	
	private void allocate(int size) {
		
		x = new double[size];
		y = new double[size];
		w = new double[size];
		
		references = new int[size];
		
		links = new int[size][];
		
	}
	
	protected HOPLayer(HOPLayer parent) {
		
		this.parent = parent;
		
	}
	
	public int size() {
		return references.length;
	}
	
	public HOPPart get(int part) {
		
		if (part < 0 || part >= size())
			return null;
		
		return new HOPPart(part);
		
	}
	
	public HOPPartIterator getIterator() {
		return new HOPPartIterator();
	}
	
}
