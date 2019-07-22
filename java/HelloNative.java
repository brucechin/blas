package blas.java;

public class HelloNative{
		static{		
			System.loadLibrary("HelloNative");
		}
		public long ptr;
		
		public static native long createCNative();
		public static native long freeCNative(long ptr);
		public static native long sayHello();
		public HelloNative(){
			System.out.println("construct");
				ptr = createCNative();
			}
		public void finalize(){
			System.out.println("deconstruct");	
			freeCNative(ptr);
			}
			
		public static void main(String[]args){
     			{
	        
				new HelloNative().sayHello();
				
				}
				System.gc();
	        }
			


}
