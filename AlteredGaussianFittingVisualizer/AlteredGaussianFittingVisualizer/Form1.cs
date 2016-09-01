using System.Runtime.InteropServices.WindowsRuntime;
using System.Windows.Forms;
using System.Windows.Forms.DataVisualization.Charting;
using DotNetMatrix;
using Net.Kniaz.LMA;
using TestLMA;

namespace AlteredGaussianFittingVisualizer
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
            string[] header;
            int startIndex;
            int endIndex;
            double[][] data = AlteredGaussian.readSkylineCopyData(@"C:\Users\mhorowitzgelb\Documents\ComplexInterference.txt",
                29.5, 31.3, out header, out startIndex, out endIndex);


            var function = new AlteredGaussian.SummedAlteredGaussian(3,header.Length,10,data);

            var lma = new LMA(function, new double[]
            {
                10, 10, 50, 10, 100, 10, 29.85, 0.01,0.15,
                20, 3000, 30, 500, 200, 300, 30.1, 0.01, 0.15,
                300, 400, 350, 50,  50, 50, 30.7, 0.01, 0.4


            }, function.AlteredData, null,new GeneralMatrix(27,27),0,300);

            lma.Fit();

            var sums = new double[data[0].Length];
            for (int i = 0; i < header.Length; i++)
            {
                for (int j = 0; j < data[0].Length; j++)
                {
                    sums[j] += data[1 + i][j];
                    
                }
                
            }

            Series dataSeries = new Series("total");
            dataSeries.ChartType = SeriesChartType.Line;

            for (int i = 0; i < sums.Length; i++)
            {
                dataSeries.Points.Add(new DataPoint(data[0][i], sums[i]));
            }

            Series fitSeries = new Series("total(fitted)");
            fitSeries.ChartType = SeriesChartType.Line;

            for (double t = 29.5; t < 31.31; t += 0.01)
            {
                var sum = 0.0;
                for (int i = 0; i < data.Length -1; i++)
                {
                    sum += function.GetY(t, i, lma.Parameters);
                }
                var point = new DataPoint(t, sum);
                fitSeries.Points.Add(point);
            }
            dataSeries.BorderWidth = 5;
            fitSeries.BorderWidth = 5;
            this.chart1.Series.Add(dataSeries);
            this.chart1.Series.Add(fitSeries);
         
        }
    }
}
