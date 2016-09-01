
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using GraphVizWrapper;
using GraphVizWrapper.Commands;
using GraphVizWrapper.Queries;

namespace TestCode
{
    class Program
    {
        static void Main(string[] args)
        {
            var rts = new Dictionary<string, ReportRow[]>();
            int peptideIndex = -1;
            int proteinIndex = -1;
            int replicateIndex = -1;
            int bestRTIndex = -1;
            int rtIndex = -1;
            var chargeIndex = -1;
            var minIndex = -1;
            var maxIndex = -1;


            using (var manReader = new StreamReader(@"C:\Users\mhorowitzgelb\Desktop\newManual.csv"))
            {
                var headerStr = manReader.ReadLine();
                string[] header = headerStr.Split(',');
                for (int i = 0; i < header.Length; i++)
                {
                    string val = header[i];
                    if (val.Contains("Peptide Modified Sequence"))
                    {
                        peptideIndex = i;
                    }
                    else if (val.Contains("Protein Name"))
                    {
                        proteinIndex = i;
                    }
                    else if (val.Contains("Replicate Name"))
                    {
                        replicateIndex = i;
                    }
                    else if (val.Equals("Retention Time"))
                    {
                        rtIndex = i;
                    }
                    else if (val.Contains("Precursor Charge"))
                    {
                        chargeIndex = i;
                    }
                    else if (val.Contains("Best Retention Time"))
                    {
                        bestRTIndex = i;
                    }
                    else if (val.Equals("Min Retention Time"))
                    {
                        minIndex = i;
                    }
                    else if (val.Equals("Max Retention Time"))
                    {
                        maxIndex = i;
                    }
                }   
                while (!manReader.EndOfStream)
                {
                    string[] line = manReader.ReadLine().Split(',');
                    var seq = line[peptideIndex];
                    var prot = line[proteinIndex];
                    var rep = line[replicateIndex];
                    var charge = line[chargeIndex];
                    var minString = line[minIndex];
                    var maxString = line[maxIndex];


                    var key = seq + prot + rep + charge;
                    
                    var bestRTString = line[bestRTIndex];
                    var rtString = line[rtIndex];
                    var bestRT = -1.0;
                    var rt = -1.0;
                    var minRt = -1.0;
                    var maxRT = -1.0;
                    if (!rtString.Contains("N/A"))
                    {
                        bestRT = Double.Parse(bestRTString);
                        minRt = Double.Parse(minString);
                        maxRT = Double.Parse(maxString);
                    }
                    if (!rtString.Contains("N/A"))
                    {
                       rt = Double.Parse(rtString);
                    }
                    
                    
                    if (rts.ContainsKey(key))
                    {
                        rts[key][0].indRTs.Add(rt);
                    }
                    else
                    {
                        var row = new ReportRow();
                        row.indRTs = new List<double>();
                        row.indRTs.Add(rt);
                        row.peptide = seq;
                        row.protein = prot;
                        row.replicate = rep;
                        row.charge = charge;
                        row.retentionTime = bestRT;
                        row.minRT = minRt;
                        row.maxRT = maxRT;
                        rts.Add(key, new[] {row, null});
                    }
                }
            }


            var doubleValues = 0;
             minIndex = -1;
             maxIndex = -1;
            int indRTIndex = -1;
            using (var tricReader = new StreamReader(@"C:\Users\mhorowitzgelb\Desktop\TRICLibraryCorrelation.csv"))
            {
                string headerStr = tricReader.ReadLine();
                string[] header = headerStr.Split(',');
                for (int i = 0; i < header.Length; i++)
                {
                    string val = header[i];
                    if (val.Contains("Peptide Modified Sequence"))
                    {
                        peptideIndex = i;
                    }
                    else if (val.Contains("Protein Name"))
                    {
                        proteinIndex = i;
                    }
                    else if (val.Contains("Replicate Name"))
                    {
                        replicateIndex = i;
                    }
                    else if (val.Equals("Best Retention Time"))
                    {
                        rtIndex = i;
                    }
                    else if (val.Contains("Precursor Charge"))
                    {
                        chargeIndex = i;
                    }
                    else if (val.Contains("Min Start Time"))
                    {
                        minIndex = i;
                    }
                    else if (val.Contains("Max End Time"))
                    {
                        maxIndex = i;
                    }
                    else if (val.Equals("Retention Time"))
                    {
                        indRTIndex = i;
                    }
                }

                while (!tricReader.EndOfStream)
                {
                    string[] line = tricReader.ReadLine().Split(',');
                    var seq = line[peptideIndex];
                    var prot = line[proteinIndex];
                    var rep = line[replicateIndex];
                    var charge = line[chargeIndex];
                    var rtString = line[rtIndex];
                    var minString = line[minIndex];
                    var maxString = line[maxIndex];
                    var inRTString = line[indRTIndex];
                    var minTime = 0.0;
                    var maxTime = 0.0;
                    var rt = 0.0;
                    var indRT = 0.0;
                    if (rtString.Contains("N/A"))
                    {
                        rt = -1.0;
                        minTime = -1.0;
                        maxTime = -1.0;
                    }
                    else
                    {
                        rt = Double.Parse(rtString);
                        minTime = Double.Parse(minString);
                        maxTime = Double.Parse(maxString);
                        indRT = Double.Parse(inRTString);
                    }
                    var hashString = seq + prot + rep+charge;
                    
                    if (!rts.ContainsKey(hashString))
                    {
                        var row = new ReportRow();
                        row.protein = prot;
                        row.peptide = seq;
                        row.charge = charge;
                        row.retentionTime = rt;
                        row.replicate = rep;
                        row.minRT = minTime;
                        row.maxRT = maxTime;
                        row.indRTs = new List<double>();
                        row.indRTs.Add(indRT);
                        rts.Add(hashString, new[] {null, row});
                    }
                    else
                    {
                        var row = rts[hashString][1];
                        if ( row != null)
                        {
                            doubleValues ++;
                            row.doubleValue = true;
                            row.indRTs.Add(indRT);
                        }
                        else {
                            row = new ReportRow();
                            row.protein = prot;
                            row.peptide = seq;
                            row.charge = charge;
                            row.retentionTime = rt;
                            row.replicate = rep;
                            row.minRT = minTime;
                            row.maxRT = maxTime;
                            row.indRTs = new List<double>();
                            row.indRTs.Add(indRT);
                            rts[hashString][1] = row;
                        }
                    }
                } 
            }

            var nonExistentTransitions = 0;
            var correctTransitions = 0;
            var outOfRangeTransitions = 0;
            var missedTransitions = 0;
            var notInTric = 0;


            var missed = new List<ReportRow[]>();
            var correct = new List<ReportRow[]>();
            var outOfRange = new List<ReportRow[]>();
            var doubleValueList = new List<ReportRow[]>();

            foreach (var pair in rts.Values)
            {
                if (pair[0] != null && pair[1] != null)
                {
                    if (pair[0].minRT == -1 && pair[1].minRT == -1)
                    {
                        correctTransitions ++;
                        correct.Add(pair);
                    }
                    else if (/*(pair[0].minRT >= pair[1].minRT && pair[0].minRT < pair[1].maxRT) || (pair[0].maxRT > pair[1].minRT && pair[0].maxRT <= pair[1].maxRT)
                        || (pair[1].minRT >= pair[0].minRT && pair[1].maxRT<= pair[0].maxRT)*/
                        Math.Abs(pair[0].retentionTime - pair[1].retentionTime) < 0.34
                        || Math.Abs(pair[0].median() - pair[1].median()) < 0.34
                        || overlap(pair[1].minRT,pair[1].maxRT, pair[0].minRT, pair[0].maxRT) >= 0.8
                        )
                    {
                        correctTransitions ++;
                        correct.Add(pair);
                    }
                    else if (pair[0].retentionTime == -1.0)
                    {
                        nonExistentTransitions ++;
                    }
                    else if (pair[1].retentionTime == -1.0)
                    {
                        missedTransitions ++;
                        missed.Add(pair);
                    }
                    else
                    {
                        outOfRangeTransitions ++;
                        outOfRange.Add(pair);
                    }
                }
                else if (pair[0]!= null)
                {
                    //missedTransitions ++;
                    //missed.Add(pair);
                    notInTric ++;
                }
                else if (pair[1] != null)
                {
                    if (pair[1].doubleValue)
                    {
                        doubleValueList.Add(pair);
                    }
                    nonExistentTransitions ++;
                }
                else
                {
                    throw new Exception("This should not happen!");
                }
            }

           
            foreach (var error in outOfRange)
            {
                var manual = error[0];
                var tric = error[1];
                Console.WriteLine(manual.retentionTime+tric.retentionTime);
            }
            foreach (var miss in missed)
            {
                var manual = miss[0];
                var tric = miss[1];
                if (tric == null)
                    continue;
                Console.WriteLine(manual.retentionTime+tric.retentionTime);
            }
        }

        static double overlap(double minA, double maxA, double minB, double maxB)
        {
            double bWidth = maxB - maxB;
            return Math.Max(0, Math.Min(maxA, maxB) - Math.Max(minA, minB))/bWidth;
        }
    }

    public class ReportRow
    {
        public double retentionTime;
        public string protein;
        public string peptide;
        public string charge;
        public string replicate;
        public bool doubleValue = false;
        public double minRT;
        public double maxRT;
        public List<double> indRTs;
        public double? cacheMedian;

        public double median()
        {
            if (cacheMedian.HasValue)
            {
                return cacheMedian.Value;
            }
            if (retentionTime == -1.0)
            {
                cacheMedian = -1.0;
                return cacheMedian.Value;
            }
            indRTs.RemoveAll(g => g == -1.0);
            if (indRTs.Count == 0)
            {
                cacheMedian = -1.0;
                return cacheMedian.Value;
            }
            cacheMedian = indRTs.OrderByDescending(g => g).ElementAt(indRTs.Count/2);
            return cacheMedian.Value;
        }
    }
}