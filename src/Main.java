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
