import ij.process.ImageProcessor;
//Na podstawie http://rosettacode.org/wiki/Image_convolution

public class Convolution {

	private static int bound(int value, int endIndex)
	  {
	    if (value < 0)
	      return 0;
	    if (value < endIndex)
	      return value;
	    return endIndex - 1;
	  }
	 
	  public static float[][] convolute(ImageProcessor inputData, ImageProcessor kernel)
	  {
	    int inputWidth = inputData.getWidth();
	    int inputHeight = inputData.getHeight();
	    int kernelWidth = kernel.getWidth();
	    int kernelHeight = kernel.getHeight();
	    if ((kernelWidth <= 0) || ((kernelWidth & 1) != 1))
	      throw new IllegalArgumentException("Kernel must have odd width");
	    if ((kernelHeight <= 0) || ((kernelHeight & 1) != 1))
	      throw new IllegalArgumentException("Kernel must have odd height");
	    int kernelWidthRadius = kernelWidth >>> 1;
	    int kernelHeightRadius = kernelHeight >>> 1;
	 
	    float[][] outputData = new float[inputWidth][inputHeight];
	    for (int i = inputWidth - 1; i >= 0; i--)
	    {
	      for (int j = inputHeight - 1; j >= 0; j--)
	      {
	        double newValue = 0.0;
	        for (int kw = kernelWidth - 1; kw >= 0; kw--){
	          for (int kh = kernelHeight - 1; kh >= 0; kh--){
	            newValue += kernel.get(kw, kh) * inputData.get(bound(i + kw - kernelWidthRadius, inputWidth),
	                          bound(j + kh - kernelHeightRadius, inputHeight));
	            outputData[i][j] = (float) newValue;
	          }
	        }
	      }
	    }
	    return outputData;
	  }
}
