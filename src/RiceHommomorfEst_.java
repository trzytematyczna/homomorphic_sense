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
		ImageProcessor Sigma_n2 = lpf(Sigma_n, (float)lpfSNR);
		
		M1 = filter2b(createImage(5, 5, 1.0/25.0), In);
		
		if(SNR == null){
			SNR = divide(M2, Sigma_n);
		}
		
		ImageProcessor LocalMean;
		ImageProcessor Rn ;
		ImageProcessor LPF1;
		ImageProcessor LPF2;
		ImageProcessor Fc1;

		double psi = -0.5772;
		double exp_psi_div2 = 1.334568251529384;

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
				LocalMean = createImage(In.getWidth(), In.getHeight(), 0.0);
			}
			
			Rn = absdiff(In, LocalMean);
				
			ImageProcessor lRn = add(multiply(Rn, compare(Rn, 0.0, NEQ)), multiply(compare(Rn, 0.0, EQ),0.001));
			lRn.log();
			LPF2 = lpf(lRn, (float)LPF);
			System.out.println("correct_rice_gauss...");
			Fc1 = correct_rice_gauss(SNR);
			LPF1 = substract(LPF2,Fc1);
			LPF1 = lpf(LPF1, (float)lpfRice, 2);
			LPF1.exp();
			ImageProcessor Mapa1 = LPF1;

//			MapaR=Mapa1.*2./sqrt(2).*exp(-psi(1)./2);
			MapaR = multiply(divide(multiply(Mapa1, 2.0), Math.sqrt(2.0)), exp_psi_div2);
//		}
		
			//GAUSSIAN!!
//		else if (noiseType == GAUSSIAN){
			System.out.println("GAUSSIAN...");
			Rn = absdiff(In, M1);
//			lRn=log(Rn.*(Rn~=0)+0.001.*(Rn==0));
			lRn = add(multiply(Rn, compare(Rn, 0.0, NEQ)), multiply(compare(Rn, 0.0, EQ),0.001));
			lRn.log();
			LPF2 = lpf(lRn,(float)LPF);
			LPF2.exp();
			ImageProcessor Mapa2 = LPF2;
//			MapaG=Mapa2.*2./sqrt(2).*exp(-psi(1)./2);
			MapaG = multiply(divide(multiply(Mapa2, 2.0), Math.sqrt(2.0)), exp_psi_div2);
			
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
			ImageProcessor h = fspecial(I.getWidth() * 2, I.getHeight() , sigma);
			
			if (Mx == 1 ||My == 1) {
				ImageProcessor lRnF = fftshift(fft(I));
				ImageProcessor lRnF2 = multiply(lRnF, h);
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
				ImageProcessor lRnF = dct(I);
				ImageProcessor lRnF2 = multiply(lRnF, h); //TODO complex multipy dimentions
				ImageProcessor If = idct(lRnF2);
				return If;
			}
			else {
				ImageProcessor lRnF = dct2(I);
				ImageProcessor lRnF2 = multiply(lRnF, h);
				ImageProcessor If = idct2(lRnF2);
				return If;
			}
		}
		return null;
	}

	public static ImageProcessor approxI1_I0(ImageProcessor z) {
//		cont=(z<1.5);
		ImageProcessor cont = compare(z, 1.5, LT);
//		z8=8.*z;
		ImageProcessor z8 = multiply(z, 8.0);

//		Mn=1-3./z8-15./2./(z8).^2-(3*5*21)./6./(z8).^3;
		ImageProcessor Mn = substract(substract(substract(1.0, divide(3.0, z8)), divide(15.0/2.0, pow(z8, 2.0))), divide((3.0*5.0*21.0)/6.0, pow(z8, 3.0)));
		
//		Md=1+1./z8+9./2./(z8).^2+(25*9)./6./(z8).^3;
		ImageProcessor Md = add(add(add(divide(1.0, z8), 1.0), divide(9.0/2.0, pow(z8,2.0))), divide((25.0*9.0)/6.0, pow(z8, 3.0)));
		
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
		K = find(z,0.0, EQ);
		
		M = applyfromfind(M, K,0.0);
		
		return M;
	}

	public static ImageProcessor besseli(int besval, ImageProcessor mat) {
		float[][] outf = new float[mat.getWidth()][mat.getHeight()]; 
		if(besval == 0){
			for(int i=0; i<mat.getWidth();i++){
				for(int j=0; j<mat.getHeight();j++){
					outf[i][j] = (float)modBesselFirstZero(mat.getPixelValue(i, j));
				}
			}
		}
		if(besval == 1){
			for(int i=0; i<mat.getWidth();i++){
				for(int j=0; j<mat.getHeight();j++){
					outf[i][j] = (float)modBesselFirstOne(mat.getPixelValue(i, j));
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
		ImageProcessor out = new FloatProcessor(outf);
		return out;
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
		            case 1:  if (mat.getPixelValue(i, j) < (float) value){//		int LT =1;
		            			temp[iterator++]=linear_index;
		            		 }
		            		 break;
		            case 2:  if (mat.getPixelValue(i, j) > (float) value){//		int	GT= 2;
		            			temp[iterator++]=linear_index;
		            		 }
		            		 break;
		            case 3:  if (mat.getPixelValue(i, j) == (float) value){//		int	EQ =3;
		            			temp[iterator++]=linear_index;
		            		 }
		            		 break;
		            case 4:  if (mat.getPixelValue(i, j) <= (float) value){//		int	LOE = 4;
		            			temp[iterator++]=linear_index;
		            		 }
		            		 break;
		            case 5:  if (mat.getPixelValue(i, j) >= (float) value){//		int	GOE =5;
		            			temp[iterator++]=linear_index;
		            		 }
		            		 break;
					case 6:  if (mat.getPixelValue(i, j) != (float) value){//		int	NEQ =6;
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
		ImageProcessor a_k = sqrt(sqrt(max(substract(multiply(pow(filter2b(Mask, pow(In, 2)), 2), 2.0), filter2b(Mask, pow(In, 4))), 0.0)));
		
//		sigma_k2=0.5.*max(filter2B(Mask,In.^2)-a_k.^2,0.01);
		ImageProcessor sigma_k2 = multiply(max(substract(filter2b(Mask, pow(In, 2)), pow(a_k, 2)), 0.01), 0.5);
		
		for(int i=1; i<N; i++){
			
//			a_k=max(filter2B(Mask,approxI1_I0(a_k.*In./sigma_k2).*In),0);
			a_k = max(filter2b(Mask, multiply(approxI1_I0(divide(multiply(a_k, In), sigma_k2)), In)), 0.0);
			
//			sigma_k2=max(0.5.*filter2B(Mask,abs(In).^2)-a_k.^2./2,0.01);
			sigma_k2 = max(substract(multiply(filter2b(Mask, pow(abs(In), 2)), 0.5), divide(pow(a_k, 2), 2.0)), 0.01);
			
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
				kernelMatrix[y][x] = (float) Math.exp(-(((x-x0)*(x-x0))/(2.0*sigma*sigma) + ((y-y0)*(y-y0))/(2.0*sigma*sigma)));
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
			pixelsAndZeros[i] = (float) 0.0;
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
			pixelsAndZeros[i] = (float) 0.0;
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
		float max = Float.MIN_VALUE;
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
		float min = Float.MAX_VALUE;
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
		matout.multiply(1.0/value);
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
		            case 1:  if (mat.getPixelValue(i, j) < (float) value){//		int LT =1;
		            			outf[i][j] = 1.0f;
		            		 }
		            		 else {
		         				outf[i][j] = 0.0f;
		            		 }
		            		 break;
		            case 2:  if (mat.getPixelValue(i, j) > (float) value){//		int	GT= 2;
		            			outf[i][j] = 1.0f;
		            		 }
		            		 else {
		         				outf[i][j] = 0.0f;
		            		 }
		            		 break;
		            case 3:  if (mat.getPixelValue(i, j) == (float) value){//		int	EQ =3;
		            			outf[i][j] = 1.0f;
		            		 }
		            		 else {
		         				outf[i][j] = 0.0f;
		            		 }
		            		 break;
		            case 4:  if (mat.getPixelValue(i, j) <= (float) value){//		int	LOE = 4;
		            			outf[i][j] = 1.0f;
		            		 }
		            		 else {
		         				outf[i][j] = 0.0f;
		            		 }
		            		 break;
		            case 5:  if (mat.getPixelValue(i, j) >= (float) value){//		int	GOE =5;
		            			outf[i][j] = 1.0f;
		            		 }
		            		 else {
		         				outf[i][j] = 0.0f;
		            		 }
		            		 break;
					case 6:  if (mat.getPixelValue(i, j) != (float) value){//		int	NEQ =6;
				    			outf[i][j] = 1.0f;
				    		 }
				    		 else {
				 				outf[i][j] = 0.0f;
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

	/**
	* The special function math library.
	* This class cannot be subclassed or instantiated because all methods are static.
	* @version 1.0
	* @author Mark Hale
	* http://web.mit.edu/fluids-modules/www/potential_flows/WavesDefBU/JSci/maths/SpecialMath.java
	*/

	
    /**
    * Modified Bessel function of first kind, order zero.
    * Based on the NETLIB Fortran function besi0 written by W. Fullerton.
    */
    public static double modBesselFirstZero(double x) {
            double y=Math.abs(x);
            if(y>3.0)
                    return Math.exp(y)*expModBesselFirstZero(x);
            else
                    return 2.75+chebyshev(y*y/4.5-1.0,bi0cs);
    }
    /**
    * Exponential scaled modified Bessel function of first kind, order zero.
    * Based on the NETLIB Fortran function besi0e written by W. Fullerton.
    */
    private static double expModBesselFirstZero(double x) {
            final double y=Math.abs(x);
            if(y>3.0) {
                    if(y>8.0)
                            return (0.375+chebyshev(16.0/y-1.0,ai02cs))/Math.sqrt(y);
                    else
                            return (0.375+chebyshev((48.0/y-11.0)/5.0,ai0cs))/Math.sqrt(y);
            } else
                    return Math.exp(-y)*(2.75+chebyshev(y*y/4.5-1.0,bi0cs));
    }
    /**
     * Modified Bessel function of first kind, order one.
     * Based on the NETLIB Fortran function besi1 written by W. Fullerton.
     */
     public static double modBesselFirstOne(double x) {
             final double y=Math.abs(x);
             if(y>3.0)
                     return Math.exp(y)*expModBesselFirstOne(x);
             else if(y==0.0)
                     return 0.0;
             else
                     return x*(0.875+chebyshev(y*y/4.5-1.0,bi1cs));
     }
     /**
     * Exponential scaled modified Bessel function of first kind, order one.
     * Based on the NETLIB Fortran function besi1e written by W. Fullerton.
     */
     private static double expModBesselFirstOne(double x) {
             final double y=Math.abs(x);
             if(y>3.0) {
                     if(y>8.0)
                             return x/y*(0.375+chebyshev(16.0/y-1.0,ai12cs))/Math.sqrt(y);
                     else
                             return x/y*(0.375+chebyshev((48.0/y-11.0)/5.0,ai1cs))/Math.sqrt(y);
             } else if(y==0.0)
                     return 0.0;
             else
                     return Math.exp(-y)*x*(0.875+chebyshev(y*y/4.5-1.0,bi1cs));
     }
    /**
     * Evaluates a Chebyshev series.
     * @param x value at which to evaluate series
     * @param series the coefficients of the series
     */
     public static double chebyshev(double x, double series[]) {
             double twox,b0=0.0,b1=0.0,b2=0.0;
             twox=2*x;
             for(int i=series.length-1;i>-1;i--) {
                     b2=b1;
                     b1=b0;
                     b0=twox*b1-b2+series[i];
             }
             return 0.5*(b0-b2);
     }
    
    private final static double bi0cs[]={
        -0.07660547252839144951,
        1.927337953993808270,
        0.2282644586920301339,
        0.01304891466707290428,
        0.00043442709008164874,
        0.00000942265768600193,
        0.00000014340062895106,
        0.00000000161384906966,
        0.00000000001396650044,
        0.00000000000009579451,
        0.00000000000000053339,
        0.00000000000000000245};
    
    private final static double ai02cs[]={
        0.05449041101410882,
        0.00336911647825569,
        0.00006889758346918,
        0.00000289137052082,
        0.00000020489185893,
        0.00000002266668991,
        0.00000000339623203,
        0.00000000049406022,
        0.00000000001188914,
        -0.00000000003149915,
        -0.00000000001321580,
        -0.00000000000179419,
        0.00000000000071801,
        0.00000000000038529,
        0.00000000000001539,
        -0.00000000000004151,
        -0.00000000000000954,
        0.00000000000000382,
        0.00000000000000176,
        -0.00000000000000034,
        -0.00000000000000027,
        0.00000000000000003};
    
    private final static double ai0cs[]={
        0.07575994494023796,
        0.00759138081082334,
        0.00041531313389237,
        0.00001070076463439,
        -0.00000790117997921,
        -0.00000078261435014,
        0.00000027838499429,
        0.00000000825247260,
        -0.00000001204463945,
        0.00000000155964859,
        0.00000000022925563,
        -0.00000000011916228,
        0.00000000001757854,
        0.00000000000112822,
        -0.00000000000114684,
        0.00000000000027155,
        -0.00000000000002415,
        -0.00000000000000608,
        0.00000000000000314,
        -0.00000000000000071,
        0.00000000000000007};
    
    private final static double bi1cs[]={
        -0.001971713261099859,
        0.40734887667546481,
        0.034838994299959456,
        0.001545394556300123,
        0.000041888521098377,
        0.000000764902676483,
        0.000000010042493924,
        0.000000000099322077,
        0.000000000000766380,
        0.000000000000004741,
        0.000000000000000024};
    
    private final static double ai12cs[]={
        0.02857623501828014,
        -0.00976109749136147,
        -0.00011058893876263,
        -0.00000388256480887,
        -0.00000025122362377,
        -0.00000002631468847,
        -0.00000000383538039,
        -0.00000000055897433,
        -0.00000000001897495,
        0.00000000003252602,
        0.00000000001412580,
        0.00000000000203564,
        -0.00000000000071985,
        -0.00000000000040836,
        -0.00000000000002101,
        0.00000000000004273,
        0.00000000000001041,
        -0.00000000000000382,
        -0.00000000000000186,
        0.00000000000000033,
        0.00000000000000028,
        -0.00000000000000003};
    
    private final static double ai1cs[]={
        -0.02846744181881479,
        -0.01922953231443221,
        -0.00061151858579437,
        -0.00002069971253350,
        0.00000858561914581,
        0.00000104949824671,
        -0.00000029183389184,
        -0.00000001559378146,
        0.00000001318012367,
        -0.00000000144842341,
        -0.00000000029085122,
        0.00000000012663889,
        -0.00000000001664947,
        -0.00000000000166665,
        0.00000000000124260,
        -0.00000000000027315,
        0.00000000000002023,
        0.00000000000000730,
        -0.00000000000000333,
        0.00000000000000071,
        -0.00000000000000006};
    
}
