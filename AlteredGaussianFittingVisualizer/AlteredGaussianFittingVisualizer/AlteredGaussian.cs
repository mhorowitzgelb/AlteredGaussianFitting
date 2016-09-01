using System;
using System.Collections.Generic;
using System.IO;
using Net.Kniaz.LMA;

namespace TestLMA
{
    class AlteredGaussian
    {
        static void run(string[] args)
        {
                        

            /*
            var t0 = 103.6374;
            using (var writer = new StreamWriter(@"C:\Users\mhorowitzgelb\Desktop\basicAsymetricGaussian.csv"))
            {
                for (double x = -3; x <= 3.1; x += 0.01)
                {
                    var _x = x + t0;
                    var y = function.GetY(_x, new[] {56564, t0, 0.00001, 1});
                    writer.WriteLine(_x+","+y);
                }
            }*/


            
            var unalteredData = loadCsv(@"C:\Users\mhorowitzgelb\Desktop\multiplePeaks.csv", 6);
            unalteredData = new[]
            {unalteredData[0], unalteredData[1]};
           
            var function = new SummedAlteredGaussian(3,1,10,unalteredData);
            //var function = new MultiTransitionAlteredGaussian(5,10,unalteredData);
            //{14000,5,5,5,5,5,43.9,0.01,0.15,1000,1300,2100,1900,900,500,44.9,0.01,0.15,100,5,4000,5,5,5,46,0.01,0.15}
            var lma = new LMA(function, new []{14000,43.9,0.01,0.15,1000,44.9,0.01,0.15,100,46,0.01,0.15},function.AlteredData,null, new DotNetMatrix.GeneralMatrix(12,12),1,1000);
            lma.Fit();
            writeCSV(lma,function,unalteredData,@"C:\Users\mhorowitzgelb\Desktop\multiPeaksFit.csv", 1);
        }

        public static double[][] readSkylineCopyData(String path, double minTime, double maxTime, out String[] header, out int startIndex, out int endIndex)
        {
            double[][] ans;
            using(var reader = new StreamReader(path))
            {
                string headerLine = reader.ReadLine();
                string[] headerSplit = headerLine.Split('\t');
                List<string> headerStrings = new List<string>();
                foreach (var s in headerSplit)
                {
                    if(s.Length > 0)
                        headerStrings.Add(s);
                }
                header = headerStrings.ToArray();
                reader.ReadLine();
                int transitions = headerStrings.Count;
                List<double>[] arr = new List<double>[transitions + 1];
                for (int i = 0; i < arr.Length; i++)
                {
                    arr[i] = new List<double>();
                }
                string line;
                startIndex = -1;
                int curIndex = -1;
                endIndex = -1;
                while ((line = reader.ReadLine()) != null)
                {
                    curIndex++;
                    string[] split = line.Split('\t');
                    if (split.Length != transitions + 1)
                        break;
                    
                    double[] vals = new double[transitions +1];
                    for (int i = 0; i < transitions + 1; i++)
                    {
                        if (!double.TryParse(split[i], out vals[i]))
                        {
                            endIndex = curIndex-1;
                            goto End;
                        }
                    }
                    if(vals[0] < minTime)
                        continue;
                    if (startIndex == -1)
                        startIndex = curIndex;


                    if (vals[0] > maxTime)
                    {
                        endIndex = curIndex;
                        break;
                    }
                    for (int i = 0; i < vals.Length; i++)
                    {
                        arr[i].Add(vals[i]);
                    }
                }
                End:
                ans = new double[transitions+1][];
                for (int i = 0; i < transitions + 1; i++)
                {
                    ans[i] = arr[i].ToArray();
                }
            }
            return ans;
        }

        public static void writeCSV(LMA lma, SummedAlteredGaussian function, double[][] originalData,string path, int transitions)
        {
            using (var writer = new StreamWriter(path))
            {
                var xList = originalData[0];
                for (int i = 0; i < xList.Length; i ++)
                {
                    var x = xList[i];
                    writer.Write(""+x);
                    for (int j = 0; j < transitions; j++)
                    {
                        var y = originalData[1 + j][i];
                        writer.Write(","+y);
                    }
                    for (int j = 0; j < transitions; j++)
                    {
                        var fitValue = function.GetY(x, j, lma.Parameters);
                        writer.Write(","+fitValue);
                    }
                    writer.WriteLine();
                }
            }
        }
        public static double[][] loadCsv(string path, int transitions)
        {
            List<double> xList = new List<double>();
            var yList = new List<List<double>>();
            for(int i = 0; i < transitions; i ++)
                yList.Add(new List<double>());
            using (var reader = new StreamReader(path))
            {
                while (!reader.EndOfStream)
                {
                    string line = reader.ReadLine();
                    var split = line.Split(',');
                    var x = double.Parse(split[0]);
                    xList.Add(x);
                    for (int i = 1; i <= transitions; i++)
                    {
                        var y = double.Parse(split[i]);
                        yList[i-1].Add(y);
                    }
                }
            }
            double[][] data = new double[1+transitions][];
            data[0] = xList.ToArray();
            for (int i = 1; i <= transitions; i++)
            {
                data[i] = yList[i - 1].ToArray();
            }
            return data;
        }


        public class SummedAlteredGaussian : LMAFunction 
        {
            private int _numGaussians;
            private int _numTransitions;
            private MultiTransitionAlteredGaussian _gaussian;
            private int _parametersPerGaussian;
            public SummedAlteredGaussian(int numGaussians, int numTransitions, double offsetLength, double[][] unalteredData)
            {
                _numGaussians = numGaussians;
                _numTransitions = numTransitions;
                _gaussian = new MultiTransitionAlteredGaussian(numTransitions, offsetLength, unalteredData);
                _parametersPerGaussian = numTransitions + 3;
            }

            public double[][] AlteredData { get { return _gaussian.AlteredData; } }

            public override double GetY(double x, double[] a)
            {
                double sum = 0;
                for (int i = 0; i < _numGaussians; i ++){
                    double[] gaussianParams = new double[_parametersPerGaussian];
                    for (int j = 0; j < _parametersPerGaussian; j ++)
                    {
                        gaussianParams[j] = a[_parametersPerGaussian*i + j];
                    }
                    sum += _gaussian.GetY(x, gaussianParams);
                }
                return sum;
            }

            public override double GetPartialDerivative(double x, double[] a, int parameterIndex)
            {
                double sum = 0;
                for (int i = 0; i < _numGaussians; i++)
                {
                    double[] gaussianParams = new double[_parametersPerGaussian];
                    for (int j = 0; j < _parametersPerGaussian; j++)
                    {
                        gaussianParams[j] = a[_parametersPerGaussian * i + j];
                    }
                    sum += _gaussian.GetPartialDerivative(x, gaussianParams, parameterIndex%_parametersPerGaussian);
                }
                return sum;
            }

            public double GetY(double x, int transition, double[] a)
            {
                double sum = 0;
                for (int i = 0; i < _numGaussians; i++)
                {
                    double[] gaussianParams = new double[_parametersPerGaussian];
                    for (int j = 0; j < _parametersPerGaussian; j++)
                    {
                        gaussianParams[j] = a[_parametersPerGaussian * i + j];
                    }
                    sum += _gaussian.GetY(x, transition, gaussianParams);
                }
                return sum;
            }
        }

        public class MultiTransitionAlteredGaussian : LMAFunction
        {
            public double[][] AlteredData { get; private set; }

            private double _start;
            private AlteredGaussianFunction _alteredGaussianFunction;
            private double _offsetLength;
            private int _transitions;
            
            public MultiTransitionAlteredGaussian(int transitions, double offsetLength, double[][] unAlterdData)
            {
                AlteredData = getAlteredData(unAlterdData, transitions, offsetLength);
                _start = unAlterdData[0][0];
                _alteredGaussianFunction = new AlteredGaussianFunction();
                _offsetLength = offsetLength;
                _transitions = transitions;
            }

            private double[][] getAlteredData(double[][] unAlterdData, int transitions, double offsetLength)
            {
                int curIndex = 0;
                var alteredData = new double[2][];
                alteredData[0] = new double[transitions* unAlterdData[0].Length];
                alteredData[1] = new double[transitions* unAlterdData[0].Length];
                for (int i = 0; i < transitions; i++)
                {
                    for (int j = 0; j < unAlterdData[0].Length; j++)
                    {
                        alteredData[0][curIndex] = unAlterdData[0][j] + i*offsetLength;
                        alteredData[1][curIndex] = unAlterdData[1 + i][j];
                        curIndex ++;
                    }
                }
                return alteredData;
            }

            public double GetY(double x, int transition, double[] a)
            {
                return _alteredGaussianFunction.GetY(x,
                    new[] {a[transition], a[_transitions], a[_transitions + 1], a[_transitions + 2]});
            }

            public override double GetY(double x, double[] a)
            {
                double? unAlteredX;
                double[] singleGaussianParams;
                GetUnAlteredValues(x,a,out unAlteredX, out singleGaussianParams);
                if (!unAlteredX.HasValue)
                    return 0;
                return _alteredGaussianFunction.GetY(unAlteredX.Value, singleGaussianParams);
            }

            

            public override double GetPartialDerivative(double x, double[] a, int parameterIndex)
            {
                double? unAlteredX;
                double[] singleGaussianParams;
                GetUnAlteredValues(x, a, out unAlteredX, out singleGaussianParams);
                if (!unAlteredX.HasValue)
                    return 0;
                if (parameterIndex < _transitions)
                {
                    return _alteredGaussianFunction.GetPartialDerivative(unAlteredX.Value, singleGaussianParams, 0);
                }
                else
                {
                    return _alteredGaussianFunction.GetPartialDerivative(unAlteredX.Value, singleGaussianParams,
                        parameterIndex - _transitions + 1);
                }
            }

            private void GetUnAlteredValues(double x, double[] a, out double? unalteredX, out double[] singleGaussianParams)
            {
                int aIndex = (int)((x - _start) / _offsetLength);
                if (aIndex < 0 || aIndex >= _transitions)
                {
                    unalteredX = null;
                    singleGaussianParams = null;
                    return;
                }
                    
                double A = a[aIndex];
                double t0 = a[_transitions];
                double mu = a[_transitions + 1];
                double sigma = a[_transitions + 2];
                unalteredX = _start + (x - _start) % _offsetLength;
                singleGaussianParams = new[] {A, t0, mu, sigma};
            }
        }

        class AlteredGaussianFunction : LMAFunction
        {
            public const int A_INDEX = 0;
            public const int T0_INDEX = 1;
            public const int MU_INDEX = 2;
            public const int SIGMA_INDEX = 3;
            public override double GetY(double t, double[] param)
            {
                var A = param[A_INDEX];
                var t0 = param[T0_INDEX];
                var mu = param[MU_INDEX];
                var sigma = param[SIGMA_INDEX];
                return ((Math.Sign(mu)*t)>=(Math.Sign(mu))*(t0-sigma*sigma/mu)) ? 
                    ((Math.Exp(-(t - t0)*(t - t0)/(2*(sigma*sigma + mu*(t - t0))))*Math.Abs(A*sigma))/
                    Math.Sqrt(Math.Abs(sigma*sigma + mu*(t - t0)))) : 0;
            }

            public override double GetPartialDerivative(double t, double[] param, int parameterIndex)
            {
                var A = param[A_INDEX];
                var t0 = param[T0_INDEX];
                var mu = param[MU_INDEX];
                var sigma = param[SIGMA_INDEX];

                if (Math.Sign(mu) * t < Math.Sign(mu) * (t0 - sigma*sigma/mu))
                {
                    return 0;
                }
                switch (parameterIndex)
                {
                    case A_INDEX:
                        return (sigma*Math.Exp(-(t - t0)*(t-t0)/(2*(sigma*sigma + mu*(t - t0))))*Math.Sign(A*sigma))/
                               Math.Sqrt(Math.Abs(sigma*sigma + mu*(t - t0)));
                    case T0_INDEX:
                        return (Math.Exp(-(t - t0) *(t-t0)/(2*(sigma*sigma + mu*(t - t0))))*Math.Abs(A*sigma)*
                                ((2*t - 2*t0)/(2*(sigma *sigma + mu*(t - t0))) -
                                 (mu * (t - t0) * (t - t0)) / (2 * (sigma * sigma + mu * (t - t0)) * (sigma * sigma + mu * (t - t0))))) / Math.Sqrt(Math.Abs(sigma *sigma + mu * (t - t0)))
                               +
                               (mu*Math.Sign(sigma *sigma + mu*(t - t0))*Math.Exp(-(t - t0) *(t-t0)/(2*(sigma *sigma + mu*(t - t0))))*
                                Math.Abs(A*sigma))/(2*Math.Pow(Math.Abs(sigma *sigma + mu*(t - t0)),3.0/2));
                    case MU_INDEX:
                        return (Math.Exp(-(t - t0)*(t - t0)/(2*(sigma*sigma + mu*(t - t0))))*Math.Abs(A*sigma)*(t - t0)*
                                (t - t0)*(t - t0))/
                               (2*Math.Sqrt(Math.Abs(sigma*sigma + mu*(t - t0)))*(sigma*sigma + mu*(t - t0))*
                                (sigma*sigma + mu*(t - t0))) -
                               (Math.Sign(sigma*sigma + mu*(t - t0))*
                                Math.Exp(-(t - t0)*(t - t0)/(2*(sigma*sigma + mu*(t - t0))))*Math.Abs(A*sigma)*(t - t0))/
                               (2*Math.Pow(Math.Abs(sigma*sigma + mu*(t - t0)), (1.5)));
                    case SIGMA_INDEX:
                        return (A*Math.Exp(-(t - t0)*(t - t0)/(2*(sigma*sigma + mu*(t - t0))))*Math.Sign(A*sigma))/
                               Math.Sqrt(Math.Abs(sigma*sigma + mu*(t - t0))) -
                               (sigma*Math.Sign(sigma*sigma + mu*(t - t0))*
                                Math.Exp(-(t - t0)*(t - t0)/(2*(sigma*sigma + mu*(t - t0))))*Math.Abs(A*sigma))/
                               Math.Pow(Math.Abs(sigma*sigma + mu*(t - t0)), (1.5)) +
                               (sigma*Math.Exp(-(t - t0)*(t - t0)/(2*(sigma*sigma + mu*(t - t0))))*Math.Abs(A*sigma)*
                                (t - t0)*(t - t0))/
                               (Math.Sqrt(Math.Abs(sigma *sigma + mu*(t - t0)))*(sigma*sigma + mu*(t - t0))*
                                (sigma*sigma + mu*(t - t0)));
                    default:
                        return double.NaN;
                }
            }

        }
    }
}
