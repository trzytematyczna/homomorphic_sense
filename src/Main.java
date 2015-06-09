import edu.emory.mathcs.jtransforms.dct.FloatDCT_1D;
import edu.emory.mathcs.jtransforms.dct.FloatDCT_2D;
import edu.emory.mathcs.jtransforms.fft.FloatFFT_1D;
import edu.emory.mathcs.jtransforms.fft.FloatFFT_2D;
import ij.ImagePlus;
import ij.plugin.GaussianBlur3D;
import ij.plugin.TextReader;
import ij.plugin.filter.GaussianBlur;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;


public class Main {
	
	public static void main(String[] args) {
		
//		System.loadLibrary("opencv_java300");
		
		TextReader textReader = new TextReader();
		ImageProcessor mriIp = textReader.open("res/MR_noisy.csv");
		ImageProcessor snrIp = textReader.open("res/MR_SNR.csv");
//		
		
//		GaussianBlur gb = new GaussianBlur();
//		float[][] kernel = gb.makeGaussianKernel(3.4, 0.0001, 256);
//		for (int i = 0; i < kernel[0].length; i++) {
//			System.out.printf("%.2f", kernel[0][i]);
//			System.out.println(" ");
//		}
		
		ImageProcessor ip = fspecial(256,256,3.4f);
		System.out.println(ip.getPixelValue(128, 128));
		System.out.println(ip.getPixelValue(50, 50));
		System.out.println(ip.getPixelValue(100, 100));
		System.out.println(ip.getPixelValue(50, 100));
		System.out.println(ip.getPixelValue(100, 50));
		
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
//	   ImageProcessor asdasd = new FloatProcessor(ppp);
//	   System.out.println(asdasd.getPixelValue(0, 0));
//	   System.out.println(asdasd.getHeight());
//	   System.out.println(asdasd.getWidth());
//		Mat mri = new Mat(mriIp.getHeight(), mriIp.getWidth(), CvType.CV_64FC1);
//		Mat snr = new Mat(snrIp.getHeight(), snrIp.getWidth(), CvType.CV_64FC1);
//		for (int i = 0; i < mriIp.getHeight(); i++){
//			for (int j = 0; j < mriIp.getWidth(); j++) {
//				mri.put(i, j, mriIp.getPixelValue(i, j));
//				snr.put(i, j, snrIp.getPixelValue(i, j));
//			}
//		}
		
//		Mat ricianMap =  RiceHommomorfEst.rice_hommomorf_est(mri, snr, 4.8, 2, RiceHommomorfEst.RICIAN, 3);
//		
//		Imgcodecs.imwrite("MR_Rician_map.csv", ricianMap);
		
		
		
	}
	
	private static ImageProcessor fspecial(int width, int height, float sigma) {
		float[][] kernelMatrix = new float[height][width];
		int x0 = width / 2;
		int y0 = width / 2;
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				kernelMatrix[y][x] = (float) Math.exp(-(((x-x0)*(x-x0))/(2*sigma*sigma) + ((y-y0)*(y-y0))/(2*sigma*sigma)));
			}
		}
		int kernelLength = height * width;
		float[] kernel = new float[kernelLength];
		int i = 0;
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < height; x++) {
				kernel[i++] = kernelMatrix[y][x];
			}
		}
		ImageProcessor kernelIP = new FloatProcessor(width, height, kernel);
		return kernelIP;
	}
	
}
