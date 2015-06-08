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
		
		for (int i = 0; i < 5; i++) {
			System.out.println(((float[])mriIp.getPixels())[i]);
		}
		System.out.println(mriIp.getPixelValue(0, 1));	
		
		
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
		
		ImageProcessor ii = fft(mriIp);
		ii = fftshift(ii);
		ii = fftshift(ii);
		ImageProcessor ff = ifft(ii);
		
		new ImagePlus("fft",ff).show();
		
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
	
	private static ImageProcessor fft(ImageProcessor ip) {
		float[] pixelsAndZeros = new float[ip.getHeight() * ip.getWidth() * 2];
		float[] pixels = (float[])ip.getPixels();
		for (int i = 0; i < ip.getHeight() * ip.getWidth(); i++) {
			pixelsAndZeros[i] = pixels[i];
		}
		for (int i = ip.getHeight() * ip.getWidth(); i < ip.getHeight() * ip.getWidth() * 2; i++) {
			pixelsAndZeros[i] = 0;
		}
		FloatFFT_2D fft2d = new FloatFFT_2D(ip.getHeight(), ip.getWidth());
		fft2d.realForwardFull(pixelsAndZeros);
		ImageProcessor fftout = new FloatProcessor(ip.getWidth() * 2, ip.getHeight(), pixelsAndZeros);
		return fftout;
	}
	
	private static ImageProcessor ifft(ImageProcessor ip) {
		float[] pixelsCopy = (float[])ip.getPixelsCopy();
		FloatFFT_2D fft2d = new FloatFFT_2D(ip.getHeight(), ip.getWidth() / 2);
		fft2d.complexInverse(pixelsCopy, true);
		float[] pixels = new float[ip.getHeight() * ip.getWidth() / 2];
		for (int i = 0; i < ip.getHeight() * ip.getWidth() / 2; i++) {
			pixels[i] = pixelsCopy[i * 2];
		}
		ImageProcessor ifftout = new FloatProcessor(ip.getWidth() / 2, ip.getHeight(), pixels);
		return ifftout;
	}
	
	private static ImageProcessor fftshift(ImageProcessor mat) {
		
		float[][] shift = new float[mat.getWidth()][mat.getHeight()];
		float[][] pixelsCopy = mat.getFloatArray();
		
		for (int i = 0; i < mat.getWidth() / 2; i++) {
			for (int j = 0; j < mat.getHeight() / 2; j++) {
				shift[i][j] = pixelsCopy[i + mat.getWidth() / 2][j + mat.getHeight() / 2];
				shift[i + mat.getWidth() / 2][j] = pixelsCopy[i][j + mat.getHeight() / 2];
				shift[i][j + mat.getHeight() / 2] = pixelsCopy[i + mat.getWidth() / 2][j];
				shift[i + mat.getWidth() / 2][j + mat.getHeight() / 2] = pixelsCopy[i][j];
			}
		}
		
		ImageProcessor ipshift = new FloatProcessor(shift);
		return ipshift;
		
	}
	
}
