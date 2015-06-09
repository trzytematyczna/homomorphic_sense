import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Properties;

import sun.security.krb5.internal.PAData;
import edu.emory.mathcs.jtransforms.dct.FloatDCT_1D;
import edu.emory.mathcs.jtransforms.dct.FloatDCT_2D;
import ij.plugin.GaussianBlur3D;
import ij.plugin.TextReader;
import ij.plugin.filter.GaussianBlur;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;


public class Main {
	public static void main(String[] args) {
		
		RiceHommomorfEst_ rhe = new RiceHommomorfEst_();
		rhe.setup("homm", null);
		rhe.run(null);
		
//		float[] a = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
//		ImageProcessor I = new FloatProcessor(3, 4, a);
//		I = RiceHommomorfEst_.padarray(I, 2, 2);
//		System.out.println(I.getWidth() + " " + I.getHeight());
//		for (int i = 0; i < I.getWidth(); i++) {
//			for (int j = 0; j < I.getHeight(); j++) {
//				System.out.print(I.getPixelValue(i, j) + " ");
//			}
//			System.out.println();
//		}
		
//		for (int i = 0; i < 5; i++) {
//			System.out.println(((float[])mriIp.getPixels())[i]);
//		}
//		System.out.println(mriIp.getPixelValue(0, 1));	
		
//		GaussianBlur gb = new GaussianBlur();
//		float[][] kernel = gb.makeGaussianKernel(3.4, 0.0001, 256);
//		for (int i = 0; i < kernel[0].length; i++) {
//			System.out.printf("%.2f", kernel[0][i]);
//			System.out.println(" ");
//		}
		
//		ImageProcessor ip = RiceHommomorfEst_.fspecial(256,256,3.4f);
//		System.out.println(ip.getPixelValue(128, 128));
//		System.out.println(ip.getPixelValue(50, 50));
//		System.out.println(ip.getPixelValue(100, 100));
//		System.out.println(ip.getPixelValue(50, 100));
//		System.out.println(ip.getPixelValue(100, 50));
		
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
		
//		float[] a = {1, 2, 3, 4, 5, 6, 7, 8};
//		ImageProcessor tst = new FloatProcessor(8, 1, a);
//		
//		ImageProcessor ii = dct(tst);
//		
//		for (int i = 0; i < ii.getPixelCount(); i++) {
//			System.out.printf("%.2f", ((float[])ii.getPixels())[i]);
//			System.out.print(" ");
//		} 
//		
//		System.out.println();
//		
//		ImageProcessor ff = idct(ii);
//
//		for (int i = 0; i < ff.getPixelCount(); i++) {
//			System.out.printf("%.2f", ((float[])ff.getPixels())[i]);
//			System.out.print(" ");
//		}
		
//		float[] a = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 2, 3, 4, 5};
//		ImageProcessor tst = new FloatProcessor(4, 4, a);
//		float[] b = {1, 2, 3, 4, 5, 6 ,7, 8, 9};
//		multiply(tst.convolve(b, 3, 3), sum(b));
//		
//		
//		for (int i = 0; i < tst.getPixelCount(); i++) {
//			System.out.printf("%.2f", ((float[])tst.getPixels())[i]);
//			System.out.print(" ");
//		}
		
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
	
}
