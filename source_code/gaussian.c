
/************************Blur(smooth) Image using gaussian 3*3 kernel ****************************************************

Function Name :  SmoothImage
Input         :  int *** , int , int , int
Output        :  int ****

*****************************************************************************************/

int*** SmoothImage(int ***XImage ,int width, int height, int channel)
{

	
     double sigma =  1.85; //(n/2 -1)*0.3 +0.8  { n = 9 ,no. of elements}
	 double conv[3][3];
	 double hg = 0;

	//Convolution kernel
	for(int i=0; i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			
			int u=i-1; //subtract from centre element index in 3*3 (1,1)
			int v=j-1;
			
			conv[i][j] =exp ( - ((u*u) + (v*v)) / (2*sigma*sigma) ); 
			hg += conv[i][j];
			
                      
		}
	}
	for(int i=0; i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			
		
			conv[i][j] = conv[i][j] / hg;
			
                      
		}
	}


     int*** sXImage = 0;
     sXImage =  CreateImageMatrix( sXImage , width , height , channel );// allocating a 3d array

	//Assigning weights
		
	 for(int i =0; i < height ; i++)
		{
			for(int j =0; j < width ; j++)
			{
			
				for(int k =0; k < 3 ; k++)
				{
					double val = 0;
					double valw =0;
					
								if(j-1 > 0  && i-1 > 0)
								{
									val += conv[0][0] * XImage[i-1][j-1][k];
									valw += conv[0][0];
								}
								if(i-1 > 0 )
								{
									val += conv[0][1] * XImage[i-1][j][k];
									valw += conv[0][1];

								}
								if(i-1 > 0 && j+1 < width)
								{
									val += conv[0][2] * XImage[i-1][j+1][k];
									valw += conv[0][2];
								}

								if(j-1 > 0 )
								{
									val += conv[1][0] * XImage[i][j-1][k];
									valw += conv[1][0];
								}

								val += conv[1][1] * XImage[i][j][k];
								valw += conv[1][1];

								if(j+1 < width)
								{
									val += conv[1][2] * XImage[i][j+1][k];
									valw += conv[1][2];
								}

								if(j-1 > 0  && i+1 < height)
								{
									val += conv[2][0] * XImage[i+1][j-1][k];
									valw += conv[2][0];
								}
								if(i+1 < height)
								{
									val += conv[2][1] * XImage[i+1][j][k];
									valw += conv[2][1];
								}
								if(j+1 < width  && i+1 < height)
								{
									val += conv[2][2] * XImage[i+1][j+1][k];
									valw += conv[2][2];
								}
								
						
							
					   	sXImage[i][j][k] = val / valw;

					
						
				}
			}
	 }


			return( sXImage);
		
		
}