import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;

import ij.*;
import ij.gui.GenericDialog;
import ij.io.FileSaver;
import ij.io.OpenDialog;
import ij.measure.ResultsTable;
import ij.plugin.TextReader;
import ij.plugin.filter.PlugInFilter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import edu.emory.mathcs.jtransforms.dct.FloatDCT_1D;
import edu.emory.mathcs.jtransforms.dct.FloatDCT_2D;
import edu.emory.mathcs.jtransforms.fft.*;

import java.util.Properties;

public class RiceHommomorfEst_ implements PlugInFilter {

	static int ex_filter_type; //# 1 - local mean, # 2 - expectation-maximization (EM).
	static int ex_window_size; //# window size for E{X},
	static int ex_iterations;  //   # number of iterations of the EM algorithm (used only by EM),
	static double lpf_f;    //        # sigma for LPF filter,
	static double lpf_f_SNR;  //      # sigma for LPF filter; used to smooth sigma(x) in SNR,
	static double lpf_f_Rice;       //# sigma for LPF filter; used to smooth Rician corrected noise map,
	static String input_filename; //# Noisy MR image,
	static String input_filenameSNR; //                  # Noisy MR image,
	static String output_filename_Gaussian; //# estimated noise map for Gaussian case,
	static String output_filename_Rician; //    # estimated noise map for Rician case,
	
	//operatory porownan
	static int LT =1; //<
	static int	GT= 2; //>
	static int	EQ =3; //==
	static int	LOE = 4; // <=
	static int	GOE =5; // >=
	static int	NEQ =6; //!=
	
	@Override
	public void run(ImageProcessor arg0) {
		
		System.out.println("Running...");
		
		System.out.println("Reading images...");
		TextReader textReader = new TextReader();
		ImageProcessor mriIp = textReader.open(input_filename);
		ImageProcessor snrIp;
		if(input_filenameSNR.equalsIgnoreCase("0")){
			snrIp = null;			
		}
		else{
			snrIp = textReader.open(input_filenameSNR);
		}
		
		System.out.println("rice_hommomorf_est...");
		long start = System.currentTimeMillis();
		ImageProcessor[] res = rice_hommomorf_est(mriIp, snrIp, lpf_f,lpf_f_SNR, lpf_f_Rice, 
					ex_filter_type,ex_window_size, ex_iterations);
		long finish = System.currentTimeMillis();
		
		String text = "Execution time: "+(finish-start)+" ms";
		GenericDialog gd = new GenericDialog("New Image");
		gd.addMessage(text);
		gd.showDialog();

		ResultsTable rician =  ResultsTable.createTableFromImage(res[0]);
		ResultsTable gauss =  ResultsTable.createTableFromImage(res[1]);
		
		System.out.println("Saving images...");
		ImagePlus rician_plus = new ImagePlus("rician", res[0]);
		FileSaver fs_rician = new FileSaver(rician_plus);
		ImagePlus gauss_plus = new ImagePlus("gauss", res[1]);
		FileSaver fs_gauss = new FileSaver(gauss_plus);
		
		fs_rician.saveAsText(new File(output_filename_Rician).getAbsolutePath());
		fs_gauss.saveAsText(new File(output_filename_Gaussian).getAbsolutePath());
		
		
		System.out.println("Showing images...");
		rician_plus.show();
		gauss_plus.show();
		
//		try {
//			rician.saveAs(output_filename_Rician);
//			gauss.saveAs(output_filename_Gaussian);
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
		
	}

	@Override
	public int setup(String arg0, ImagePlus arg1) {
		try {
			System.out.print("Running setup...");
			
   		    OpenDialog dialog = new OpenDialog("Choose properties file!");
		    String path_prop =dialog.getPath();
			File file = new File(path_prop);
			FileInputStream fileInput = new FileInputStream(file);
			Properties properties = new Properties();
			properties.load(fileInput);
			fileInput.close();
			ex_filter_type = Integer.parseInt(properties.getProperty("ex_filter_type"));
			ex_window_size = Integer.parseInt(properties.getProperty("ex_window_size"));
			ex_iterations = Integer.parseInt(properties.getProperty("ex_iterations"));
			lpf_f = Double.parseDouble(properties.getProperty("lpf_f"));
			lpf_f_SNR = Double.parseDouble(properties.getProperty("lpf_f_SNR"));
			lpf_f_Rice = Double.parseDouble(properties.getProperty("lpf_f_Rice"));
			input_filename = properties.getProperty("input_filename");
			input_filename = input_filename.replace("'", "");
			input_filenameSNR = properties.getProperty("input_filenameSNR");
			input_filenameSNR = input_filenameSNR.replace("'", "");
			output_filename_Gaussian = properties.getProperty("output_filename_Gaussian");
			output_filename_Gaussian = output_filename_Gaussian.replace("'", "");
			output_filename_Rician = properties.getProperty("output_filename_Rician");
			output_filename_Rician = output_filename_Rician.replace("'", "");
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		System.out.println("DONE");
		return  DOES_ALL+NO_IMAGE_REQUIRED;
	}
	

//	lpfSNR dla lpf(sigma_n, ...)
//	lpf dla wywolywania gauss + rician
// 	lpfRice do uzycia przy lpf() po correct_rice_gauss
	
	public static ImageProcessor[] rice_hommomorf_est(ImageProcessor In, ImageProcessor SNR, 
			double LPF, double lpfSNR, double lpfRice, int Modo,int winsize, int iteration_no){
		
		ImageProcessor[] result = new FloatProcessor[2];
		
		ImageProcessor MapaR = null;
		ImageProcessor MapaG = null;
		
		int[] Ws = {winsize, winsize};
		
		System.out.println("em_ml_rice2D...");
		ImageProcessor [] em_ml = em_ml_rice2D(In, iteration_no, Ws);
		ImageProcessor M2 = em_ml[0];
		ImageProcessor Sigma_n = em_ml[1];
		ImageProcessor M1;
		Sigma_n = lpf(Sigma_n, (float)lpfSNR);
		
		M1 = filter2b(createImage(5, 5, 0.25), In);
		
		if(SNR.getHeight() == 1 && SNR.getHeight() == 0){
			SNR = divide(M2, Sigma_n);
		}
		
		ImageProcessor LocalMean;
		ImageProcessor Rn ;
		ImageProcessor LPF1;
		ImageProcessor LPF2;
		ImageProcessor Fc1;

		double psi = -0.5772;
		double exp_psi_div2 = 1.3346;

		//RICIAN!!
//		if(noiseType == RICIAN){
		System.out.println("RICIAN...");
			if(Modo == 1){
				LocalMean = M1;
			}
			else if (Modo == 2){
				LocalMean = M2;
			}
			else{
				LocalMean = createImage(In.getWidth(), In.getHeight(), 0);
			}
			
			Rn = absdiff(In, LocalMean);
				
			ImageProcessor lRn = add(multiply(Rn, compare(Rn, 0, NEQ)), multiply(compare(Rn, 0, EQ),0.001));
			lRn.log();
			LPF2 = lpf(lRn, (float)LPF);
			System.out.println("correct_rice_gauss...");
			Fc1 = correct_rice_gauss(SNR);
			LPF1 = substract(LPF2,Fc1);
			LPF1 = lpf(LPF1, (float)lpfRice, 2);
			LPF1.exp();
			ImageProcessor Mapa1 = LPF1;

			MapaR = multiply(Mapa1, 2);
			MapaR.sqrt();
			MapaR = multiply(MapaR, exp_psi_div2);
//		}
		
			//GAUSSIAN!!
//		else if (noiseType == GAUSSIAN){
			System.out.println("GAUSSIAN...");
			Rn = absdiff(In, M1);
			lRn = add(multiply(Rn, compare(Rn, 0, NEQ)), multiply(compare(Rn, 0, EQ),0.001));
			lRn.log();
			LPF2 = lpf(lRn,(float)LPF);
			LPF2.exp();
			ImageProcessor Mapa2 = LPF2;
			MapaG = multiply(Mapa2, 2);
			MapaG.sqrt();
			MapaG = multiply(MapaG, exp_psi_div2);
			
//		}
		result[0] = MapaR;
		result[1] = MapaG;
		
		System.out.println("DONE");
		return result;
	}

	public static ImageProcessor correct_rice_gauss(ImageProcessor SNR) {
		//Fc=Coefs(1)+Coefs(2).*a1+Coefs(3).*a1.^2+Coefs(4).*a1.^3+Coefs(5).*a1.^4+Coefs(6).*a1.^5+Coefs(7).*a1.^6+Coefs(8).*a1.^7+Coefs(9).*a1.^8;
		//Fc=Fc.*(a1<=7);       

		double[] Coefs = {-0.289549906258443,	-0.0388922575606330,	0.409867108141953,	-0.355237628488567,
				0.149328280945610,	-0.0357861117942093,	0.00497952893859122,	-0.000374756374477592,	1.18020229140092e-05};

		ImageProcessor Fc = add(multiply(pow(SNR, 8), Coefs[8]), multiply(pow(SNR, 7), Coefs[7]));
		Fc = add(Fc, add(multiply(pow(SNR, 6), Coefs[6]), multiply(pow(SNR, 5), Coefs[5])));
		Fc = add(Fc, add(multiply(pow(SNR, 4), Coefs[4]), multiply(pow(SNR, 3), Coefs[3])));
		Fc = add(Fc, add(multiply(pow(SNR, 2), Coefs[2]), multiply(pow(SNR, 1), Coefs[1])));
		Fc = add(Fc, Coefs[0]);
		
		Fc = multiply(Fc, compare(SNR, 7, LOE));

		return Fc;
	}

	public static ImageProcessor filter2(ImageProcessor h, ImageProcessor I) {
		int Mx = h.getHeight();
		int My = h.getWidth();
		
		if(Mx % 2 == 0 || My % 2 == 0){
			System.err.println("filter2b size of h must be odd");
			return null;
		}
		
		int Nx = (Mx - 1) / 2;
		int Ny = (My - 1) / 2;
		ImageProcessor It = padarray(I, Nx, Ny);
		ImageProcessor I2 = filter2(h, It);
//		I_out=I2((Nx+1):end-Nx,(Ny+1):end-Ny);
		ImageProcessor I_out = submatrix(I2, Ny + 1, I2.getWidth() - Ny,  Nx + 1, I2.getHeight() - Nx);
		
		return I_out;
		
	}

	public static ImageProcessor submatrix(ImageProcessor mat, int beginW, int endW,int beginH, int endH) {
		float[][] outf = new float[endW-beginW + 1][endH-beginH + 1];
		for (int i = 0; i < endW-beginW + 1; i++) {
			for (int j = 0; j < endH-beginH + 1; j++) {
				outf[i][j]= mat.getPixelValue(i + beginW,j + beginH);
			}
		}
		ImageProcessor out = new FloatProcessor(outf);
		return out;
	}

	public static ImageProcessor padarray(ImageProcessor in, int nx, int ny) {
		int rows = in.getHeight() + ny * 2;
		int cols = in.getWidth() + nx * 2;
		float[][] outf = new float[rows][cols];
		for (int y = 0; y < rows; y++) {
			for (int x = 0; x < cols; x++) {
				int cx = -1;
				int cy = -1;
				if (x < nx) cx = 0;
				else if (x < nx + in.getWidth()) cx = x - nx;
				else cx = in.getWidth() - 1;
				if (y < ny) cy = 0;
				else if (y < ny + in.getHeight()) cy = y - ny;
				else cy = in.getHeight() - 1;
				outf[y][x] = in.getPixelValue(cx, cy);
			}
		}
		ImageProcessor output = new FloatProcessor(outf);
		return output;
	}

	public static ImageProcessor lpf(ImageProcessor I, float sigma) {
		return lpf(I, sigma, 2);
	}

	public static ImageProcessor lpf(ImageProcessor I, float sigma, int MODO) {

		if (MODO == 1) {
			int Mx = I.getHeight();
			int My = I.getWidth();
			ImageProcessor h = fspecial(I.getWidth(), I.getHeight() , sigma);
			
			if (Mx == 1 ||My == 1) {
				ImageProcessor lRnF = fftshift(fft(I));
				ImageProcessor lRnF2 = multiply(lRnF, h); //TODO complex multipy dimentions
				ImageProcessor If = ifft(fftshift(lRnF2));
				return If;
			}
			else {
				ImageProcessor lRnF = fftshift(fft2(I));
				ImageProcessor lRnF2 = multiply(lRnF, h);
				ImageProcessor If = ifft2(fftshift(lRnF2));
				return If;
			}
		}
		else if (MODO == 2) {
			int Mx = I.getHeight();
			int My = I.getWidth();
			ImageProcessor h = fspecial(My * 2, Mx * 2 , sigma * 2);
			h = submatrix(h, Mx, h.getHeight() - 1, My, h.getWidth() - 1);
			
			if (Mx == 1 ||My == 1) {
				ImageProcessor lRnF = fftshift(dct(I));
				ImageProcessor lRnF2 = multiply(lRnF, h); //TODO complex multipy dimentions
				ImageProcessor If = idct(fftshift(lRnF2));
				return If;
			}
			else {
				ImageProcessor lRnF = fftshift(dct2(I));
				ImageProcessor lRnF2 = multiply(lRnF, h);
				ImageProcessor If = idct2(fftshift(lRnF2));
				return If;
			}
		}
		return null;
	}

	public static ImageProcessor approxI1_I0(ImageProcessor z) {
//		cont=(z<1.5);
		ImageProcessor cont = compare(z, 1.5, LT);
//		z8=8.*z;
		ImageProcessor z8 = multiply(z, 8);

//		Mn=1-3./z8-15./2./(z8).^2-(3*5*21)./6./(z8).^3;
		ImageProcessor Mn = substract(substract(substract(1, divide(3, z8)), divide(15/2, pow(z8, 2))), divide((3*5*21)/6, pow(z8, 3)));
		
//		Md=1+1./z8+9./2./(z8).^2+(25*9)./6./(z8).^3;
		ImageProcessor Md = add(add(add(divide(1, z8), 1), divide(9/2, pow(z8,2))), divide((25*9)/6, pow(z8, 3)));
		
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
		ImageProcessor K = null;
		if(sum(cont)>1){
			K = find(z, 1.5, LT);			
			M = applyfromfind(M, K, divide(besseli(1, valuesfromfind(z,K)), besseli(0, valuesfromfind(z,K))));
		}
		K = find(z,0, EQ);
		
		M = applyfromfind(M, K,0.0);
		
		return M;
	}

	public static ImageProcessor besseli(int besval, ImageProcessor mat) {
		float[][] outf = new float[mat.getWidth()][mat.getHeight()]; 
		if(besval == 0){
			for(int i=0; i<mat.getWidth();i++){
				for(int j=0; j<mat.getHeight();j++){
					outf[i][j] = (float) bessj0(mat.getPixelValue(i, j));
				}
			}
		}
		if(besval == 1){
			for(int i=0; i<mat.getWidth();i++){
				for(int j=0; j<mat.getHeight();j++){
					outf[i][j] =  (float) bessj1(mat.getPixelValue(i, j));
				}
			}
		}
		return new FloatProcessor(outf);
	}

	public static ImageProcessor valuesfromfind(ImageProcessor z, ImageProcessor k) {
		float[][] outf = new float[1][k.getHeight()];
		int linear_index = 0;
		int iterator = 0;
		
		for(int i=0; i<z.getWidth();i++){
			for(int j=0; j<z.getHeight();j++){
				if(linear_index == k.getPixelValue(0, iterator)){
					outf[0][iterator] = z.getPixelValue(i, j);					
					iterator++;
				}
				linear_index++;
			}
		}		
		return new FloatProcessor(outf);
	}

	public static ImageProcessor applyfromfind(ImageProcessor mat,ImageProcessor k, ImageProcessor z) {
		float linear_index = 0;
		int iterator = 0;
		
		float[] z_linear = (float [])z.getPixels();
		
		float[][] outf= new float[mat.getWidth()][mat.getHeight()];
		for(int i=0; i<mat.getWidth();i++){
			for(int j=0; j<mat.getHeight();j++){
				if(linear_index == k.getPixelValue(0, iterator)){
					outf[i][j] = z_linear[iterator];//z.getPixelValue(1, iterator);					
					iterator++;
				}
				else{
					outf[i][j] = mat.getPixelValue(i, j);
				}
				linear_index++;
			}
		}
		return new FloatProcessor(outf);
	}
	
	public static ImageProcessor applyfromfind(ImageProcessor mat,ImageProcessor k, double val) {
		float linear_index = 0;
		int iterator = 0;
		
		float[][] outf= new float[mat.getWidth()][mat.getHeight()];
		for(int i=0; i<mat.getWidth();i++){
			for(int j=0; j<mat.getHeight();j++){
				if(linear_index == k.getPixelValue(0, iterator)){
					iterator++;
					outf[i][j] = (float) val;	
				}
				else{
					outf[i][j] = mat.getPixelValue(i, j);
				}
				linear_index++;
			}
		}
		return new FloatProcessor(outf);
	}

	public static ImageProcessor find(ImageProcessor mat, double value, int cmpop) {
		int linear_index = 0;
		int iterator = 0;
		int[] temp = new int[mat.getWidth()*mat.getHeight()];
		for(int i=0; i<mat.getWidth();i++){
			for(int j=0; j<mat.getHeight();j++){
				switch (cmpop) {
		            case 1:  if (mat.getPixelValue(i, j) < value){//		int LT =1;
		            			temp[iterator++]=linear_index;
		            		 }
		            		 break;
		            case 2:  if (mat.getPixelValue(i, j) > value){//		int	GT= 2;
		            			temp[iterator++]=linear_index;
		            		 }
		            		 break;
		            case 3:  if (mat.getPixelValue(i, j) == value){//		int	EQ =3;
		            			temp[iterator++]=linear_index;
		            		 }
		            		 break;
		            case 4:  if (mat.getPixelValue(i, j) <= value){//		int	LOE = 4;
		            			temp[iterator++]=linear_index;
		            		 }
		            		 break;
		            case 5:  if (mat.getPixelValue(i, j) >= value){//		int	GOE =5;
		            			temp[iterator++]=linear_index;
		            		 }
		            		 break;
					case 6:  if (mat.getPixelValue(i, j) != value){//		int	NEQ =6;
								temp[iterator++]=linear_index;
				    		 }
				    		 break;
				}
				
				linear_index++;		 
		   }
		}
		float[][] outf = new float[1][iterator];
		for(int i = 0; i<iterator; i++){
			outf[0][i] = temp[i];
		}
		
		return new FloatProcessor(outf);
	}

	public static double sum(ImageProcessor mat1) {
		double res = 0.0;
		for(int i=0; i<mat1.getWidth();i++){
			for(int j=0; j<mat1.getHeight();j++){
				res += mat1.getPixelValue(i, j);
			}			
		}
		return res;
	}

	public static ImageProcessor[] em_ml_rice2D(ImageProcessor In, int N, int[] Ws) {
		
		int Mx = In.getHeight();
		int My = In.getWidth();
		
		float prod = Ws[0] * Ws[1];
		ImageProcessor Mask = divide(createImage(Ws[1],Ws[0], 1.0), prod);
		
//		a_k=sqrt(sqrt(max(2.*filter2B(Mask,In.^2).^2-filter2B(Mask,In.^4),0)));
		ImageProcessor a_k = sqrt(sqrt(max(substract(multiply(pow(filter2b(Mask, pow(In, 2)), 2), 2), filter2b(Mask, pow(In, 4))), 0)));
		
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
	
	public static ImageProcessor fspecial(int width, int height, float sigma) {
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
	
	public static ImageProcessor idct2(ImageProcessor ip) {
		float[] pixelsCopy = (float[]) ip.getPixelsCopy();
		FloatDCT_2D dct2d = new FloatDCT_2D(ip.getHeight(), ip.getWidth());
		dct2d.inverse(pixelsCopy, true);
		ImageProcessor dct2out = new FloatProcessor(ip.getWidth(), ip.getHeight(), pixelsCopy);
		return dct2out;
	}

	public static ImageProcessor dct2(ImageProcessor ip) {
		float[] pixelsCopy = (float[]) ip.getPixelsCopy();
		FloatDCT_2D dct2d = new FloatDCT_2D(ip.getHeight(), ip.getWidth());
		dct2d.forward(pixelsCopy, true);
		ImageProcessor dct2out = new FloatProcessor(ip.getWidth(), ip.getHeight(), pixelsCopy);
		return dct2out;
	}

	public static ImageProcessor idct(ImageProcessor ip) {
		float[] pixelsCopy = (float[]) ip.getPixelsCopy();
		FloatDCT_1D dct1d = new FloatDCT_1D(ip.getWidth());
		dct1d.inverse(pixelsCopy, true);
		ImageProcessor dct2out = new FloatProcessor(ip.getWidth(), 1, pixelsCopy);
		return dct2out;
	}

	public static ImageProcessor dct(ImageProcessor ip) {
		float[] pixelsCopy = (float[]) ip.getPixelsCopy();
		FloatDCT_1D dct1d = new FloatDCT_1D(ip.getWidth());
		dct1d.forward(pixelsCopy, true);
		ImageProcessor dct2out = new FloatProcessor(ip.getWidth(), 1, pixelsCopy);
		return dct2out;
	}
	
	public static ImageProcessor fftshift(ImageProcessor mat) {
						
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
	
	public static ImageProcessor fft(ImageProcessor ip) {
		float[] pixelsAndZeros = new float[ip.getWidth() * 2];
		float[] pixels = (float[])ip.getPixels();
		for (int i = 0; i < ip.getWidth(); i++) {
			pixelsAndZeros[i] = pixels[i];
		}
		for (int i = ip.getWidth(); i < ip.getWidth() * 2; i++) {
			pixelsAndZeros[i] = 0;
		}
		FloatFFT_1D fft1d = new FloatFFT_1D(ip.getWidth());
		fft1d.realForwardFull(pixelsAndZeros);
		ImageProcessor fftout = new FloatProcessor(ip.getWidth() * 2, 1, pixelsAndZeros);
		return fftout;
	}
	
	public static ImageProcessor ifft(ImageProcessor ip) {
		float[] pixelsCopy = (float[])ip.getPixelsCopy();
		FloatFFT_1D fft1d = new FloatFFT_1D(ip.getWidth() / 2);
		fft1d.complexInverse(pixelsCopy, true);
		float[] pixels = new float[ip.getWidth() / 2];
		for (int i = 0; i < ip.getWidth() / 2; i++) {
			pixels[i] = pixelsCopy[i * 2];
		}
		ImageProcessor ifftout = new FloatProcessor(ip.getWidth() / 2, 1, pixels);
		return ifftout;
	}
	
	public static ImageProcessor fft2(ImageProcessor ip) {
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
	
	public static ImageProcessor ifft2(ImageProcessor ip) {
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
	
	public static ImageProcessor createImage(int width, int height, double fill) {
		float[][] data = new float[width][height];
		for (int i = 0; i < width; i++) {
			for (int j = 0; j < height; j++) {
				data[i][j] = (float) fill;
			}
		}
		ImageProcessor ip = new FloatProcessor(data);
		return ip;
	}
	
	public static ImageProcessor absdiff(ImageProcessor mat1, ImageProcessor mat2) {
		return abs(substract(mat1, mat2));
	}
	

	public static ImageProcessor filter2b(ImageProcessor mask, ImageProcessor mat) {
		ImageProcessor matout = new FloatProcessor(mat.getFloatArray());
		matout.convolve((float [])mask.getPixels(), mask.getWidth(), mask.getHeight());
		matout = multiply(matout, sum(mask));
		return matout;
	}
	

	public static ImageProcessor add(ImageProcessor mat1, ImageProcessor mat2) { //done
		float[][] outf = new float[mat1.getWidth()][mat1.getHeight()];
		for(int i=0; i<mat1.getWidth();i++){
			for(int j=0; j<mat1.getHeight();j++){
				outf[i][j] = mat1.getPixelValue(i, j)+mat2.getPixelValue(i, j);
			}			
		}
		ImageProcessor out = new FloatProcessor(outf);
		return out;
	}
	
	public static ImageProcessor add(ImageProcessor mat, double value) {//done
		ImageProcessor matout = new FloatProcessor(mat.getFloatArray());
		matout.add(value);
		return matout;
	}
	
	public static ImageProcessor abs(ImageProcessor mat) {//done
		ImageProcessor matout = new FloatProcessor(mat.getFloatArray());
		matout.abs();
		return matout;
	}
	
	public static ImageProcessor substract(ImageProcessor mat1, ImageProcessor mat2) { //done
		float[][] outf = new float[mat1.getWidth()][mat1.getHeight()];
		for(int i=0; i<mat1.getWidth();i++){
			for(int j=0; j<mat1.getHeight();j++){
				outf[i][j] = mat1.getPixelValue(i, j)-mat2.getPixelValue(i, j);
			}			
		}
		ImageProcessor out = new FloatProcessor(outf);
		return out;
	}
	
	public static ImageProcessor substract(double value, ImageProcessor mat) {//done
		ImageProcessor newip = substract(createImage(mat.getWidth(),mat.getHeight(), value), mat);
		return newip ;
	}
	
	public static ImageProcessor sqrt(ImageProcessor mat) {//done
		ImageProcessor matout = new FloatProcessor(mat.getFloatArray());
		matout.sqrt();
		return matout;
	}
	
	public static ImageProcessor pow(ImageProcessor mat, double power) {//done
		ImageProcessor matout = new FloatProcessor(mat.getFloatArray());
		ImageProcessor mattmp = new FloatProcessor(mat.getFloatArray());
		for (int i=0; i < power - 1; i++) {
			matout = multiply(matout, mattmp);			
		}
		return matout;
	}
	
	public static double max(ImageProcessor mat) { //done
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
	
	public static ImageProcessor max(ImageProcessor mat, double d) { //chyba done
		float[][] outf = new float[mat.getWidth()][mat.getHeight()];
		for(int i=0; i<mat.getWidth();i++){
			for(int j=0; j<mat.getHeight();j++){
				float matval =mat.getPixelValue(i, j);	
				if(d > matval){ //wstawiam wartosc value
					outf[i][j] = (float) d;
				}
				else{ //wstawiam ten sam element co byl
					outf[i][j] = matval;
				}
			}
		}
		ImageProcessor newip = new FloatProcessor(outf);
		return newip;
	}
	
	public static double min(ImageProcessor mat) { //done
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
	
	public static ImageProcessor divide(ImageProcessor mat, double value) {	//done
		ImageProcessor matout = new FloatProcessor(mat.getFloatArray());
		matout.multiply(1/value);
		return matout;
	}
	
	public static ImageProcessor divide(double value, ImageProcessor mat) {	//done
		ImageProcessor newip = divide(createImage(mat.getWidth(),mat.getHeight(), value), mat);
		return newip ;
	}
	
	public static ImageProcessor divide(ImageProcessor mat1, ImageProcessor mat2) {	//done
		float[][] outf = new float[mat1.getWidth()][mat1.getHeight()];
		for(int i=0; i<mat1.getWidth();i++){
			for(int j=0; j<mat1.getHeight();j++){
				outf[i][j] = mat1.getPixelValue(i, j)/mat2.getPixelValue(i, j);
			}			
		}
		ImageProcessor out = new FloatProcessor(outf);
		return out;
	}
	
	public static ImageProcessor compare(ImageProcessor mat, double value, int cmpop) { //done
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
					case 6:  if (mat.getPixelValue(i, j) != value){//		int	NEQ =6;
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
	
	public static ImageProcessor multiply(ImageProcessor mat, double value) {	//done
		ImageProcessor matout = new FloatProcessor(mat.getFloatArray());
		matout.multiply(value);
		return matout;
	}
	public static ImageProcessor multiply(ImageProcessor mat1, ImageProcessor mat2) {	//done
		float[][] outf = new float[mat1.getWidth()][mat1.getHeight()];
		for(int i=0; i<mat1.getWidth();i++){
			for(int j=0; j<mat1.getHeight();j++){
				outf[i][j] = mat1.getPixelValue(i, j)*mat2.getPixelValue(i, j);
			}			
		}
		ImageProcessor out = new FloatProcessor(outf);
		return out;
	}

//ZRODLO: http://www.atnf.csiro.au/computing/software/gipsy/sub/bessel.c
		
		static double bessj0( double x )
		/*------------------------------------------------------------*/
		/* PURPOSE: Evaluate Bessel function of first kind and order  */
		/*          0 at input x                                      */
		/*------------------------------------------------------------*/
		{
		   double ax,z;
		   double xx,y,ans,ans1,ans2;

		   if ((ax=Math.abs(x)) < 8.0) {
		      y=x*x;
		      ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7
		         +y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
		      ans2=57568490411.0+y*(1029532985.0+y*(9494680.718
		         +y*(59272.64853+y*(267.8532712+y*1.0))));
		      ans=ans1/ans2;
		   } else {
		      z=8.0/ax;
		      y=z*z;
		      xx=ax-0.785398164;
		      ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
		         +y*(-0.2073370639e-5+y*0.2093887211e-6)));
		      ans2 = -0.1562499995e-1+y*(0.1430488765e-3
		         +y*(-0.6911147651e-5+y*(0.7621095161e-6
		         -y*0.934935152e-7)));
		      ans=Math.sqrt(0.636619772/ax)*(Math.cos(xx)*ans1-z*Math.sin(xx)*ans2);
		   }
		   return ans;
		}


		static double bessj1( double x )
		/*------------------------------------------------------------*/
		/* PURPOSE: Evaluate Bessel function of first kind and order  */
		/*          1 at input x                                      */
		/*------------------------------------------------------------*/
		{
		   double ax,z;
		   double xx,y,ans,ans1,ans2;

		   if ((ax=Math.abs(x)) < 8.0) {
		      y=x*x;
		      ans1=x*(72362614232.0+y*(-7895059235.0+y*(242396853.1
		         +y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
		      ans2=144725228442.0+y*(2300535178.0+y*(18583304.74
		         +y*(99447.43394+y*(376.9991397+y*1.0))));
		      ans=ans1/ans2;
		   } else {
		      z=8.0/ax;
		      y=z*z;
		      xx=ax-2.356194491;
		      ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4
		         +y*(0.2457520174e-5+y*(-0.240337019e-6))));
		      ans2=0.04687499995+y*(-0.2002690873e-3
		         +y*(0.8449199096e-5+y*(-0.88228987e-6
		         +y*0.105787412e-6)));
		      ans=Math.sqrt(0.636619772/ax)*(Math.cos(xx)*ans1-z*Math.sin(xx)*ans2);
		      if (x < 0.0) ans = -ans;
		   }
		   return ans;
		}

}
