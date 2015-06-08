import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import ij.*;
import ij.gui.NewImage;
import ij.plugin.TextReader;
import ij.plugin.filter.PlugInFilter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.FloatProcessor;
import edu.emory.mathcs.jtransforms.fft.*;

public class RiceHommomorfEst_ implements PlugInFilter {

	static int GAUSSIAN = 1;
	static int RICIAN = 2;
	int ex_filter_type = 1;     //# 1 - local mean, 
//            # 2 - expectation-maximization (EM).
	int ex_window_size =5; //# window size for E{X},
	int ex_iterations= 10;  //   # number of iterations of the EM algorithm (used only by EM),
	double lpf_f = 3.4;    //        # sigma for LPF filter,
	double lpf_f_SNR = 1.2;  //      # sigma for LPF filter; used to smooth sigma(x) in SNR,
	double lpf_f_Rice= 5.4;       //# sigma for LPF filter; used to smooth Rician corrected noise map,
	String input_filename = "MR_noisy.csv";                    //# Noisy MR image,
	String input_filenameSNR = "MR_SNR.csv";  //                  # Noisy MR image,
	String output_filename_Gaussian = "MR_Gaussian_Map.csv";   //# estimated noise map for Gaussian case,
	String output_filename_Rician = "MR_Rician_Map.csv"  ;//    # estimated noise map for Rician case,
	
	
	int LT =1;
	int	GT= 2;
	int	EQ =3;
	int	LOE = 4;
	int	GOE =5;
	
	@Override
	public void run(ImageProcessor arg0) {
		TextReader textReader = new TextReader();
		ImageProcessor mriIp = textReader.open("res/MR_noisy.csv");
		ImageProcessor snrIp = textReader.open("res/MR_SNR.csv");
		ImagePlus asd = new ImagePlus();
	    NewImage.createFloatImage("asd", 2, 2, 2, 2);
	    float[] jedne = {1,1,1,1,1,1};
	    float[] dwa = {1,1,1,1,1,1};
	    float[][] www = {jedne,dwa};
	   ImageProcessor asdasd = new FloatProcessor(www);
//				ImagePlus("asd", mriIp);
	}

	@Override
	public int setup(String arg0, ImagePlus arg1) {
		// TODO Auto-generated method stub
		
		BufferedReader br = null;
		try {
 
			String sCurrentLine;
			br = new BufferedReader(new FileReader("C:\\testing.txt"));
			while ((sCurrentLine = br.readLine()) != null) {
				System.out.println(sCurrentLine);
			}
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				if (br != null)br.close();
			} catch (IOException ex) {
				ex.printStackTrace();
			}
		}
		
		return  DOES_ALL+NO_IMAGE_REQUIRED;
	}
	
	public static ImageProcessor rice_hommomorf_est(ImageProcessor In, ImageProcessor SNR, double LPF, int Modo, int noiseType, int winsize) {
		
		int[] Ws = {winsize, winsize};
		
		ImageProcessor [] em_ml = em_ml_rice2D(In, 10, Ws);
		ImageProcessor M2 = em_ml[0];
		ImageProcessor Sigma_n = em_ml[1];
		ImageProcessor M1;
		
		
		Sigma_n = lpf(Sigma_n,1.2);
		
		M1 = filter2b(new FloatProcessor().createProcessor(5, 5). ImageProcessor.ones(5, 5), In);
		
		if(SNR.getHeight() == 1 && SNR.getHeight() == 0){
			SNR = divide(M2, Sigma_n);
		}
		
		ImageProcessor LocalMean;
		ImageProcessor Rn;
		ImageProcessor LPF1;
		ImageProcessor LPF2;
		ImageProcessor Fc1;

		if(noiseType == RICIAN){
			if(Modo == 1){
				LocalMean = M1;
			}
			else if (Modo == 2){
				LocalMean = M2;
			}
			else{
				LocalMean = new ImageProcessor();
			}
			
			Rn = absdiff(In, M1);
					
			
//			lRn.elementLog();
//			LPF2 = lpf(lRn,LPF);
//			Fc1 = correct_rice_gauss(SNR);
//			LPF1 = LPF2.minus(Fc1);
//			LPF1 = lpf(LPF1, LPF+2,2);
//			ImageProcessor Mapa1 = LPF1.elementExp();
//			
//			double psi = -0.5772;
//			for (int i=0; i<Mapa1.numRows(); i++){
//				for (int j=0; j<Mapa1.numCols(); j++){
//					double val = Mapa1.get(i, j)*2/Math.sqrt(2)*Math.exp(-psi/2);
//					Mapa1.set(i, j, val);
//				}
//			}
			
			return null;
		}
		else if (noiseType == GAUSSIAN){
			
		}
		
		return null;
	}

	private static ImageProcessor correct_rice_gauss(ImageProcessor SNR) {

//		//Fc=Coefs(1)+Coefs(2).*a1+Coefs(3).*a1.^2+Coefs(4).*a1.^3+Coefs(5).*a1.^4+Coefs(6).*a1.^5+Coefs(7).*a1.^6+Coefs(8).*a1.^7+Coefs(9).*a1.^8;
//		//
//		//Fc=Fc.*(a1<=7);       
//
//		ImageProcessor Fc = new ImageProcessor(SNR.numRows(), SNR.numCols());
//		Fc = multipleM(SNR, coeffs[1]);
//		Fc.plus(coeffs[0]).plus(multipleM(SNR.elementPower(2), coeffs[2])).plus(multipleM(SNR.elementPower(3), coeffs[3])).plus(multipleM(SNR.elementPower(4), coeffs[4])).plus(multipleM(SNR.elementPower(5), coeffs[5])).plus(multipleM(SNR.elementPower(6), coeffs[6])).plus(multipleM(SNR.elementPower(7), coeffs[7])).plus(multipleM(SNR.elementPower(8), coeffs[8]));
//		
//		for (int i=0; i<Fc.numRows(); i++){
//			for (int j=0; j<Fc.numCols(); j++){
//				double snr_val = SNR.get(i, j);
//				if(snr_val <=7){
//					Fc.set(i, j, Fc.get(i, j)*snr_val);
//				}
//			}
//		}
		
		return null;
	}

	private static ImageProcessor filter2b(ImageProcessor h, ImageProcessor I) {
		
		int Mx = h.getHeight();
		int My = h.getWidth();
//		System.out.println(Mx);
//		System.out.println(My);
		
		if(Mx % 2 == 0 || My % 2 == 0){
			System.err.println("filter2b size of h must be odd");
			return null;
		}
		
		int Nx = (Mx - 1) / 2;
		int Ny = (My - 1) / 2;
		ImageProcessor It = padarray(I, Nx, Ny);
		
		
		ImageProcessor I2 = filter2(h, It);
//		I_out=I2((Nx+1):end-Nx,(Ny+1):end-Ny);
		ImageProcessor I_out = I2.submat(new Range(Nx + 1, I2.getHeight() - Nx), new Range(Ny + 1, I2.getWidth() - Ny));
		
		return I_out;
		
	}

	private static ImageProcessor padarray(ImageProcessor in, int nx, int ny) {
		int rows = in.getHeight() + nx * in.getHeight();
		int cols = in.getWidth() + ny * in.getWidth();
		ImageProcessor output = new ImageProcessor(rows, cols, in.type());
//		for (int i = 0; i < rows; i++) {
//			for (int j = 0; j < cols; j++) {
//				double value = in.get(i / (1 + nx), j / (1 + ny));
//				output.set(i, j, value);
//			}
//		}
		return output;
	}

	private static ImageProcessor lpf(ImageProcessor sigma_n, double d) {
		// TODO Auto-generated method stub
		return null;
	}

	private static ImageProcessor lpf(ImageProcessor I, double sigma, int MODO) {
		if (MODO == 1) {
			int Mx = I.getHeight();
			int My = I.getWidth();
			ImageProcessor h = Imgproc.getGaussianKernel(I.getHeight() , sigma, CvType.CV_64F);
			h = divide(h, max(h));
			if (Mx == 1 ||My == 1) {
				ImageProcessor lRnF = new ImageProcessor(I.size(), I.type());
				Core.dft(I, lRnF);
				fftshift(lRnF);
				ImageProcessor lRnF2 = multiply(lRnF, h);
				fftshift(lRnF2);
				ImageProcessor If = new ImageProcessor(I.size(), I.type());
				Core.idft(lRnF2, If);
				return If;
			}
			else {
				ImageProcessor lRnF = new ImageProcessor(I.size(), I.type());
				Core.dft(I, lRnF);
				fftshift(lRnF);
				ImageProcessor lRnF2 = multiply(lRnF, h);
				fftshift(lRnF2);
				ImageProcessor If = new ImageProcessor(I.size(), I.type());
				Core.idft(lRnF2, If);
				return If;
			}
		} else if (MODO == 2) {
			int Mx = I.getHeight();
			int My = I.getWidth();
			ImageProcessor h = Imgproc.getGaussianKernel(I.getHeight() * 2 , sigma * 2, CvType.CV_64F);
			h = divide(h, max(h));
			h = h.submat(new Range(Mx + 1, h.getHeight()), new Range(My + 1, h.getWidth()));
			if (Mx == 1 ||My == 1) {
				ImageProcessor lRnF = new ImageProcessor(I.size(), I.type());
				Core.dct(I, lRnF);
				ImageProcessor lRnF2 = multiply(lRnF, h);
				ImageProcessor If = new ImageProcessor(I.size(), I.type());
				Core.idct(lRnF2, If);
				return If;
			}
			else {
				ImageProcessor lRnF = new ImageProcessor(I.size(), I.type());
				Core.dct(I, lRnF);
				ImageProcessor lRnF2 = multiply(lRnF, h);
				ImageProcessor If = new ImageProcessor(I.size(), I.type());
				Core.idct(lRnF2, If);
				return If;
			}
		}
		return null;
	}
	
	private static ImageProcessor approxI1_I0(ImageProcessor z) {
//		cont=(z<1.5);
		ImageProcessor cont = compare(z, 1.5, Core.CMP_LE);
//		z8=8.*z;
		ImageProcessor z8 = multiply(z, 8);

//		Mn=1-3./z8-15./2./(z8).^2-(3*5*21)./6./(z8).^3;
		ImageProcessor Mn = substract(substract(substract(1, divide(3, z8)), divide(15/2, pow(z8, 2))), divide((3*5*21)/6, pow(z8, 3)));
		
//		Md=1+1./z8+9./2./(z8).^2+(25*9)./6./(z8).^3;
		ImageProcessor Md = add(add(add(divide(1, z8), 1), divide(9/2, z8)), divide((25*9)/6, pow(z8, 3)));
		
//		M=Mn./Md;
		ImageProcessor M = divide(Mn, Md);
//
//		%K=find(isnan(M));
//		%M(K)=besseli(1,z(K))./besseli(0,z(K));
//
//		if (sum(cont)>1)
//		K=find(z<1.5);
//		M(K)=besseli(1,z(K))./besseli(0,z(K));
//		end
//		K=find(z==0);
//		M(K)=0;
		return null;
	}

	public static ImageProcessor[] em_ml_rice2D(ImageProcessor In, int N, int[] Ws) {
		
		int Mx = In.getHeight();
		int My = In.getWidth();
		
		double prod = Ws[0] * Ws[1];
		ImageProcessor Mask = divide(createImage(Ws[1],Ws[0], 1.0), prod);
		
//		a_k=sqrt(sqrt(max(2.*filter2B(Mask,In.^2).^2-filter2B(Mask,In.^4),0)));
		ImageProcessor a_k = sqrt(max(substract(multiply(pow(filter2b(Mask, pow(In, 2)), 2), 2), filter2b(Mask, pow(In, 4))), 0));
		
//		sigma_k2=0.5.*max(filter2B(Mask,In.^2)-a_k.^2,0.01);
		ImageProcessor sigma_k2 = multiply(max(substract(filter2b(Mask, pow(In, 2)), pow(a_k, 2)), 0.01), 0.5);
		
		for(int i=1; i<N; i++){
			
//			a_k=max(filter2B(Mask,approxI1_I0(a_k.*In./sigma_k2).*In),0);
			a_k = max(filter2b(Mask, multiply(approxI1_I0(divide(multiply(a_k, In), sigma_k2)), In)), 0);
			
//			sigma_k2=max(0.5.*filter2B(Mask,abs(In).^2)-a_k.^2./2,0.01);
			sigma_k2 = max(substract(multiply(filter2b(Mask, pow(abs(In), 2)), 0.5), divide(pow(a_k, 2), 2)), 0.01);
			
		}
			
		
		ImageProcessor[] result = new ImageProcessor[2];
		result[0] = a_k;
		result[1] = sqrt(sigma_k2);
		
		return result;
	}
	
	private static void fftshift(ImageProcessor mat) {
//		ImageProcessor m1 = mat.rowRange(0, 0).colRange(mat.getHeight() / 2, mat.getWidth() / 2);
//		ImageProcessor m2 = mat.rowRange(0, mat.getWidth() / 2).colRange(mat.getHeight() / 2, mat.getWidth() / 2);
//		ImageProcessor m3 = mat.rowRange(mat.getHeight() / 2, 0).colRange(mat.getHeight() / 2, mat.getWidth() / 2);
//		ImageProcessor m4 = mat.rowRange(mat.getHeight() / 2, mat.getHeight() / 2).colRange(mat.getHeight() / 2, mat.getWidth() / 2);
		Point p1 = new Point(0, 0);
		Point p2 = new Point(0, (int)mat.getWidth() / 2);
		Point p3 = new Point((int)mat.getHeight() / 2, 0);
		Point p4 = new Point((int)mat.getHeight() / 2, (int)mat.getWidth() / 2);
		Rect r1 = new Rect(0, 0, mat.getHeight() / 2, mat.getWidth() / 2);
		Rect r2 = new Rect(0, mat.getWidth() / 2, mat.getHeight() / 2, mat.getWidth() / 2);
		Rect r3 = new Rect(mat.getHeight() / 2, 0, mat.getHeight() / 2, mat.getWidth() / 2);
		Rect r4 = new Rect(mat.getHeight() / 2, mat.getHeight() / 2, mat.getHeight() / 2, mat.getWidth() / 2);
		ImageProcessor m1 = mat.submat(r1);
		ImageProcessor m2 = mat.submat(r2); 
		ImageProcessor m3 = mat.submat(r3);
		ImageProcessor m4 = mat.submat(r4);
		
		double[] b1 = new double[(int)m1.total() * mat.channels()];
		m1.get(0, 0, b1);
		double[] b2 = new double[(int)m2.total() * mat.channels()];
		m2.get(0, 0, b2);
		double[] b3 = new double[(int)m3.total() * mat.channels()];
		m3.get(0, 0, b3);
		double[] b4 = new double[(int)m4.total() * mat.channels()];
		m4.get(0, 0, b4);
		
		mat.put((int)p1.x, (int)p1.y, b4);
		mat.put((int)p2.x, (int)p2.y, b3);
		mat.put((int)p3.x, (int)p3.y, b2);
		mat.put((int)p4.x, (int)p4.y, b1);
		
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
	
	private static ImageProcessor createImage(int width, int height, double fill) {
		float[][] data = new float[width][height];
		for (int i = 0; i < width; i++) {
			for (int j = 0; j < height; j++) {
				data[i][j] = (float) fill;
			}
		}
		ImageProcessor ip = new FloatProcessor(data);
		return ip;
	}
	
	private static ImageProcessor absdiff(ImageProcessor in, ImageProcessor m1) {
		return null;
	}
	

	private static ImageProcessor filter2(ImageProcessor mat1, ImageProcessor mat2) {
		ImageProcessor out = new ImageProcessor(mat1.getHeight(), mat1.getWidth(), mat1.type());
		Imgproc.filter2D(mat1, out, -1, mat2);
		return out;
	}
	
	
	private static ImageProcessor add(ImageProcessor mat1, ImageProcessor mat2) { //done
		float[][] outf = new float[mat1.getWidth()][mat1.getHeight()];
		for(int i=0; i<mat1.getWidth();i++){
			for(int j=0; j<mat1.getHeight();j++){
				outf[i][j] = mat1.getPixelValue(i, j)+mat2.getPixelValue(i, j);
			}			
		}
		ImageProcessor out = new FloatProcessor(outf);
		return out;
	}
	
	private static ImageProcessor add(ImageProcessor mat, double value) {//done
		mat.add(value);
		return mat;
	}
	
	private static ImageProcessor abs(ImageProcessor mat) {//done
		mat.abs();
		return mat;
	}
	
	private static ImageProcessor substract(ImageProcessor mat1, ImageProcessor mat2) { //done
		float[][] outf = new float[mat1.getWidth()][mat1.getHeight()];
		for(int i=0; i<mat1.getWidth();i++){
			for(int j=0; j<mat1.getHeight();j++){
				outf[i][j] = mat1.getPixelValue(i, j)-mat2.getPixelValue(i, j);
			}			
		}
		ImageProcessor out = new FloatProcessor(outf);
		return out;
	}
	
	private static ImageProcessor substract(double value, ImageProcessor mat) {//done
		ImageProcessor newip = substract(createImage(mat.getWidth(),mat.getHeight(), value), mat);
		return newip ;
	}
	
	private static ImageProcessor sqrt(ImageProcessor mat) {//done
		mat.sqrt();
		return mat;
	}
	
	private static ImageProcessor pow(ImageProcessor mat, double power) {//done
		for(int i=0; i<power; i++){
			mat = multiply(mat, mat);			
		}
		return mat;
	}
	
	private static double max(ImageProcessor mat) { //done
		float max = Integer.MIN_VALUE;
		for(int i=0; i<mat.getWidth();i++){
			for(int j=0; j<mat.getHeight();j++){
				float val =mat.getPixelValue(i, j);
				if(max < val){
					max = val;
				}
			}			
		}
		return max;
	}
	
	private static ImageProcessor max(ImageProcessor mat, double value) { //chyba done
		float max = Integer.MIN_VALUE;
		float[][] outf = new float[1][mat.getHeight()];
		for(int j=0; j<mat.getHeight();j++){
			for(int i=0; i<mat.getWidth();i++){
				float val =mat.getPixelValue(i, j);
				if(max < val){
					max = val;
				}
			}
			outf[1][j] = max;
			max = Integer.MIN_VALUE;
		}
		ImageProcessor newip = new FloatProcessor(outf);
		return newip;
	}
	
	private static double min(ImageProcessor mat) { //done
		float min = Integer.MAX_VALUE;
		for(int i=0; i<mat.getWidth();i++){
			for(int j=0; j<mat.getHeight();j++){
				float val =mat.getPixelValue(i, j);
				if(min > val){
					min = val;
				}
			}			
		}
		return min;
	}
	
	private static ImageProcessor divide(ImageProcessor mat, double value) {	//done
		mat.multiply(1/value);
		return mat;
	}
	
	private static ImageProcessor divide(double value, ImageProcessor mat) {	//done
		ImageProcessor newip = divide(createImage(mat.getWidth(),mat.getHeight(), value), mat);
		return newip ;
	}
	
	private static ImageProcessor divide(ImageProcessor mat1, ImageProcessor mat2) {	//done
		float[][] outf = new float[mat1.getWidth()][mat1.getHeight()];
		for(int i=0; i<mat1.getWidth();i++){
			for(int j=0; j<mat1.getHeight();j++){
				outf[i][j] = mat1.getPixelValue(i, j)/mat2.getPixelValue(i, j);
			}			
		}
		ImageProcessor out = new FloatProcessor(outf);
		return out;
	}
	
	private static ImageProcessor compare(ImageProcessor mat, double value, int cmpop) { //done
		float[][] outf = new float[mat.getWidth()][mat.getHeight()];
		for(int i=0; i<mat.getWidth();i++){
			for(int j=0; j<mat.getHeight();j++){
				switch (cmpop) {
		            case 1:  if (mat.getPixelValue(i, j) < value){//		int LT =1;
		            			outf[i][j] = 1;
		            		 }
		            		 else {
		         				outf[i][j] = 0;
		            		 }
		            		 break;
		            case 2:  if (mat.getPixelValue(i, j) > value){//		int	GT= 2;
		            			outf[i][j] = 1;
		            		 }
		            		 else {
		         				outf[i][j] = 0;
		            		 }
		            		 break;
		            case 3:  if (mat.getPixelValue(i, j) == value){//		int	EQ =3;
		            			outf[i][j] = 1;
		            		 }
		            		 else {
		         				outf[i][j] = 0;
		            		 }
		            		 break;
		            case 4:  if (mat.getPixelValue(i, j) <= value){//		int	LOE = 4;
		            			outf[i][j] = 1;
		            		 }
		            		 else {
		         				outf[i][j] = 0;
		            		 }
		            		 break;
		            case 5:  if (mat.getPixelValue(i, j) >= value){//		int	GOE =5;
		            			outf[i][j] = 1;
		            		 }
		            		 else {
		         				outf[i][j] = 0;
		            		 }
		            		 break;
				        }
			}			
		}
		ImageProcessor out = new FloatProcessor(outf);
		return out;
	}
	
	private static ImageProcessor multiply(ImageProcessor mat, double value) {	//done
		mat.multiply(value);
		return mat;
	}
	private static ImageProcessor multiply(ImageProcessor mat1, ImageProcessor mat2) {	//done
		float[][] outf = new float[mat1.getWidth()][mat1.getHeight()];
		for(int i=0; i<mat1.getWidth();i++){
			for(int j=0; j<mat1.getHeight();j++){
				outf[i][j] = mat1.getPixelValue(i, j)*mat2.getPixelValue(i, j);
			}			
		}
		ImageProcessor out = new FloatProcessor(outf);
		return out;
	}
	
}
