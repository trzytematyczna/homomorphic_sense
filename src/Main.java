import edu.emory.mathcs.jtransforms.dct.FloatDCT_1D;
import edu.emory.mathcs.jtransforms.dct.FloatDCT_2D;
import edu.emory.mathcs.jtransforms.fft.FloatFFT_1D;
import edu.emory.mathcs.jtransforms.fft.FloatFFT_2D;
import ij.ImagePlus;
import ij.plugin.TextReader;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;


public class Main {
	
	public static void main(String[] args) {
		
//		System.loadLibrary("opencv_java300");
		
		TextReader textReader = new TextReader();
		ImageProcessor mriIp = textReader.open("res/MR_noisy.csv");
		ImageProcessor snrIp = textReader.open("res/MR_SNR.csv");
//		
//		for (int i = 0; i < 5; i++) {
//			System.out.println(((float[])mriIp.getPixels())[i]);
//		}
//		System.out.println(mriIp.getPixelValue(0, 1));	
		
		
//		FloatFFT_2D fft2d = new FloatFFT_2D(4, 2);
//		fft2d.realForwardFull(a);
//		for (int i = 0; i < 16; i++) {
//			System.out.printf("%.2f", a[i]);
//			System.out.println(" ");
//		}
//		fft2d.complexInverse(a, true);
//		for (int i = 0; i < 16; i++) {
//			System.out.printf("%.2f", a[i]);
//			System.out.println(" ");
//		}
//		
//		ImageProcessor ii = fft(mriIp);
//		ii = fftshift(ii);
//		ii = fftshift(ii);
//		ImageProcessor ff = ifft(ii);
//		
//		new ImagePlus("fft",ff).show();
		
		float[] a = {1, 2, 3, 4, 5, 6, 7, 8};
		ImageProcessor tst = new FloatProcessor(8, 1, a);
		
		ImageProcessor ii = dct(tst);
		
		for (int i = 0; i < ii.getPixelCount(); i++) {
			System.out.printf("%.2f", ((float[])ii.getPixels())[i]);
			System.out.print(" ");
		} 
		
		System.out.println();
		
		ImageProcessor ff = idct(ii);

		for (int i = 0; i < ff.getPixelCount(); i++) {
			System.out.printf("%.2f", ((float[])ff.getPixels())[i]);
			System.out.print(" ");
		}
		
		
//		System.out.println(mriIp.getPixelValue(0, 0));
//		System.out.println(mriIp.getPixelValue(1, 1));
//		mriIp.add(4);
//		System.out.println(mriIp.getPixelValue(0, 0));
//		System.out.println(mriIp.getPixelValue(1, 1));
//		ImagePlus asd = new ImagePlus("ASD", mriIp);
//	    float[] jedne = {2,1,1,1,1,1};
//	    float[] dwa = {1,1,1,1,1,1};
//	    float[][] www = {jedne,dwa};
//	    float[][] ppp = new float[7][8];
//	    System.out.println(www[3][0]);
//		float[] a = {1,2,3,4};
//		float[] b = {5,6,7,8};
//		float[] c = {9,0,1,2};
//		float[] d = {3,4,5,6};
//		float[][] ppp = {a,b,c,d};
//	   ImageProcessor mat = new FloatProcessor(ppp);
//	   System.out.println(mat.getWidth());
//	   System.out.println(mat.getHeight());
//	   float[] f1={1,0,0};
//	   float[] f2={0,1,0};
//	   float[] f3={1,0,0};
//	   
//	   float[][] fil = {f1,f2,f3};
//	   
//	   ImageProcessor filter = new FloatProcessor(fil);
	   
//	   filter = RiceHommomorfEst_.rotateMatrixRight(RiceHommomorfEst_.rotateMatrixRight(filter));
//
//	   float[][] res = Convolution.convolute(mat, filter);
//	   for (int i = 0; i < res.length; i++){
//			for (int j = 0; j < res.length; j++) {
//				System.out.print(res[i][j]+",");
//			}
//			System.out.println();
//		}
//
//	   
//	   ImageProcessor filtmat = RiceHommomorfEst_.filter2(filter, mat);
//	   
//	   System.out.println("asdasdasd");
//	   float[] f = {1,0,0,0,1,0,1,0,0};
//	   
//	   mat.convolve(f,3,3);
//	   
//		for (int i = 0; i < mat.getWidth(); i++){
//			for (int j = 0; j < mat.getHeight(); j++) {
//				System.out.print(mat.getPixelValue(i, j)+",");
//			}
//			System.out.println();
//		}
//		System.out.println();
//		System.out.println();
//		for (int i = 0; i < filtmat.getWidth(); i++){
//			for (int j = 0; j < filtmat.getHeight(); j++) {
//				System.out.print(filtmat.getPixelValue(i, j)+",");
//			}
//			System.out.println();
//		}
//		System.out.println();
//		System.out.println();
//		System.out.println();
//		mat = RiceHommomorfEst_.rotateMatrixRight(RiceHommomorfEst_.rotateMatrixRight(mat));
//		for (int i = 0; i < mat.getWidth(); i++){
//			for (int j = 0; j < mat.getHeight(); j++) {
//				System.out.print(mat.getPixelValue(i, j)+",");
//			}
//			System.out.println();
//		}
		
//		Mat ricianMap =  RiceHommomorfEst.rice_hommomorf_est(mri, snr, 4.8, 2, RiceHommomorfEst.RICIAN, 3);
//		
//		Imgcodecs.imwrite("MR_Rician_map.csv", ricianMap);
		
	}
	
	
	private static ImageProcessor idct2(ImageProcessor ip) {
		float[] pixelsCopy = (float[]) ip.getPixelsCopy();
		FloatDCT_2D dct2d = new FloatDCT_2D(ip.getHeight(), ip.getWidth());
		dct2d.inverse(pixelsCopy, true);
		ImageProcessor dct2out = new FloatProcessor(ip.getWidth(), ip.getHeight(), pixelsCopy);
		return dct2out;
	}

	private static ImageProcessor dct2(ImageProcessor ip) {
		float[] pixelsCopy = (float[]) ip.getPixelsCopy();
		FloatDCT_2D dct2d = new FloatDCT_2D(ip.getHeight(), ip.getWidth());
		dct2d.forward(pixelsCopy, true);
		ImageProcessor dct2out = new FloatProcessor(ip.getWidth(), ip.getHeight(), pixelsCopy);
		return dct2out;
	}

	private static ImageProcessor idct(ImageProcessor ip) {
		float[] pixelsCopy = (float[]) ip.getPixelsCopy();
		FloatDCT_1D dct1d = new FloatDCT_1D(ip.getWidth());
		dct1d.inverse(pixelsCopy, true);
		ImageProcessor dct2out = new FloatProcessor(ip.getWidth(), 1, pixelsCopy);
		return dct2out;
	}

	private static ImageProcessor dct(ImageProcessor ip) {
		float[] pixelsCopy = (float[]) ip.getPixelsCopy();
		FloatDCT_1D dct1d = new FloatDCT_1D(ip.getWidth());
		dct1d.forward(pixelsCopy, true);
		ImageProcessor dct2out = new FloatProcessor(ip.getWidth(), 1, pixelsCopy);
		return dct2out;
	}
	
}
