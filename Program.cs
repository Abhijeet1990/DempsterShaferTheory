using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using static SwagShopLibrary.DSTheory;
using ScottPlot;

namespace TestDSTheory
{
    class Program
    {
        static void Main(string[] args)
        {
            Dictionary<string, double> hyp_val = new Dictionary<string, double> { { "ab", 0.6 }, { "bc", 0.3 }, { "a", 0.1 }, { "ad", 0.0 } };

            MassFunction m = new MassFunction(hyp_val);
            m.ComputeFocalSet();
            m.ComputeFrameOfDiscernment();
            var bel = m.belief("ab");
            var pl = m.plausibility("ab");
            var q = m.commonality("ab");
            // get belief function and then get mass function from belief function
            var bf = m.beliefFunction();
            var mf = m.fromBelief(bf);

            // get plausibility function and then get mass function from plausibility function
            var pf = m.plausFunction();
            var mf_pf = m.fromPlausibility(pf);

            // get commonality function and then get mass function from commonality function
            var qf = m.commFunction();
            var mf_qf = m.fromCommonality(qf);

            //combine deterministically
            Dictionary<string, double> hyp_val2 = new Dictionary<string, double> { { "ab", 0.4 }, { "bc", 0.4 }, { "a", 0.2 }, { "ad", 0.0 } };
            var combined = m.combinedDeterministically(hyp_val2, hyp_val, false); // conjunctive combination

            // combine disjunctive
            var disjunctive_combined = m.combinedDeterministically(hyp_val2, hyp_val, true);

            // combine with normalization
            var combined_norm = m.combine(hyp_val2, hyp_val, true, 0, false, false);

            var combinedSample = m.combinedDirectSampling(hyp_val2, hyp_val, 1000, false); // conjunctive combination
            var combinedImportanceSample = m.combinedImportanceSampling(hyp_val2, hyp_val, 1000);
            var samples = m.sample(1000, true, mf_pf);

            // test gbt
            var mf_from_gbt = m.gbt(pf, false, 0);

            var pignistic = m.pignistic();

            Dictionary<string, double> mcc1 = new Dictionary<string, double> { { "ab", 0.3 }, { "bc", 0.5 }, { "abc", 0.2 } };
            Dictionary<string, double> mcc2 = new Dictionary<string, double> { { "b", 0.3 }, { "bc", 0.4 }, { "abc", 0.3 } };
            m.combineCautious(mcc1, mcc2, false);


            //compute for 1000 samples
            int sampleSize = 1000;

            //if combined_cautious
            bool combinedCautious = false;

            //Time
            List<double> time = new List<double>();
            for (double i = 0; i < sampleSize; i++) time.Add(i + 1);

            // Now create a window for creating intrusions and causing attacks on the grid
            // Let say a cyber intrusion was performed on time [5,20] .... then [80,120]... then get control [200,250] 
            //create a list of tuple 
            List<Tuple<int, int>> cyber_intrusions_ids1 = new List<Tuple<int, int>>();
            Tuple<int, int> event1 = new Tuple<int, int>(5, 20); cyber_intrusions_ids1.Add(event1);
            Tuple<int, int> event2 = new Tuple<int, int>(80, 120); cyber_intrusions_ids1.Add(event2);
            Tuple<int, int> event3 = new Tuple<int, int>(200, 250); cyber_intrusions_ids1.Add(event3);
            Tuple<int, int> event4 = new Tuple<int, int>(400, 470); cyber_intrusions_ids1.Add(event4);
            Tuple<int, int> event5 = new Tuple<int, int>(700, 790); cyber_intrusions_ids1.Add(event5);

            List<Tuple<int, int>> cyber_intrusions_ids2 = new List<Tuple<int, int>>();
            Tuple<int, int> event6 = new Tuple<int, int>(15, 25); cyber_intrusions_ids2.Add(event6);
            Tuple<int, int> event7 = new Tuple<int, int>(90, 135); cyber_intrusions_ids2.Add(event7);
            Tuple<int, int> event8 = new Tuple<int, int>(220, 260); cyber_intrusions_ids2.Add(event8);
            Tuple<int, int> event9 = new Tuple<int, int>(390, 450); cyber_intrusions_ids2.Add(event9);
            Tuple<int, int> event10 = new Tuple<int, int>(670, 730); cyber_intrusions_ids2.Add(event10);

            List<Tuple<int, int>> cyber_intrusions_ids3 = new List<Tuple<int, int>>();
            Tuple<int, int> event11 = new Tuple<int, int>(90, 110); cyber_intrusions_ids3.Add(event11);
            Tuple<int, int> event12 = new Tuple<int, int>(120, 130); cyber_intrusions_ids3.Add(event12);
            Tuple<int, int> event13 = new Tuple<int, int>(230, 290); cyber_intrusions_ids3.Add(event13);
            Tuple<int, int> event14 = new Tuple<int, int>(380, 440); cyber_intrusions_ids3.Add(event14);
            Tuple<int, int> event15 = new Tuple<int, int>(670, 740); cyber_intrusions_ids3.Add(event15);
            // perform physical attack to modify measurements from [260,270]

            //var sensor1_mf_list = generaterandommassfunction(samplesize);
            //var sensor2_mf_list = generaterandommassfunction(samplesize);
            //var sensor3_mf_list = generaterandommassfunction(samplesize);

            var sensor1_mf_list = GenerateRandomMassFunctionScenario(sampleSize, cyber_intrusions_ids1);
            var sensor2_mf_list = GenerateRandomMassFunctionScenario(sampleSize, cyber_intrusions_ids2);
            var sensor3_mf_list = GenerateRandomMassFunctionScenario(sampleSize, cyber_intrusions_ids3);


            List<string> ranked_Hyp_ByBelief = new List<string>();
            List<double> ranked_count_ByBelief = new List<double>();

            List<string> ranked_Hyp_ByPlausibility = new List<string>();
            List<double> ranked_count_ByPlausibility = new List<double>();

            List<string> rankedPignistic = new List<string>();
            List<double> ranked_count_ByPignistic = new List<double>();

            // construct mf using Generalized bayesian Theorem and store it here for every sample fused.
            List<Dictionary<string, double>> gbt_sample_fused = new List<Dictionary<string, double>>();
            List<string> ranked_mf_gbt = new List<string>();
            List<double> ranked_count_ByGBT = new List<double>();

            List<double> conflictMeasure = new List<double>();
            List<double> hartleyMeasure = new List<double>();

            for (int i = 0; i < sampleSize; i++)
            {
                var s1_mf = sensor1_mf_list[i];
                var s2_mf = sensor2_mf_list[i];
                var fused = m.combinedDeterministically(s1_mf, s2_mf, false);
                //if (combinedCautious && (s1_mf.Count != 1) && (s2_mf.Count != 1)) fused = m.combineCautious(s1_mf, s2_mf);
                if (combinedCautious) fused = m.combineCautious(s1_mf, s2_mf,false);
                var s3_mf = sensor3_mf_list[i];
                var fused3 = m.combinedDeterministically(fused, s3_mf, false);
                //if (combinedCautious && (fused.Count != 1) && (s3_mf.Count != 1)) fused3 = m.combineCautious(fused, s3_mf);
                if (combinedCautious) fused3 = m.combineCautious(fused, s3_mf,false);

                //conflict measure
                MassFunction t1 = new MassFunction(fused);
                MassFunction t2 = new MassFunction(s3_mf);
                var conflict = m.conflict(t1, t2, 0);
                conflictMeasure.Add(conflict);

                //hartley measure
                var hartley = m.hartley_measure(fused3);
                hartleyMeasure.Add(hartley);

                MassFunction newM = new MassFunction(fused3);
                newM.ComputeFocalSet();
                newM.ComputeFrameOfDiscernment();
                var bel_func = newM.beliefFunction();
                var plaus_func = newM.plausFunction();

                // test gbt 
                var gbt_mf = newM.gbt(plaus_func, false, 0);
                gbt_sample_fused.Add(gbt_mf);
                

                // decision basis
                var max_hyp_bel = bel_func.OrderByDescending(x => x.Value).First().Key;
                ranked_Hyp_ByBelief.Add(max_hyp_bel);
                ranked_count_ByBelief.Add(Convert.ToDouble(max_hyp_bel.Length));

                var max_hyp_plaus = plaus_func.OrderByDescending(x => x.Value).First().Key;
                ranked_Hyp_ByPlausibility.Add(max_hyp_plaus);
                ranked_count_ByPlausibility.Add(Convert.ToDouble(max_hyp_plaus.Length));

                var pignisticFunction = newM.pignistic();
                if (pignisticFunction != null && pignisticFunction.Count!=0)
                {
                    var max_pignistic = pignisticFunction.OrderByDescending(x => x.Value).First().Key;
                    rankedPignistic.Add(max_pignistic);
                    ranked_count_ByPignistic.Add(Convert.ToDouble(max_pignistic.Length));
                }
                else
                {
                    rankedPignistic.Add(max_hyp_plaus);
                    ranked_count_ByPignistic.Add(Convert.ToDouble(max_hyp_plaus.Length));
                }

                var max_mf_gbt = gbt_mf.OrderByDescending(x => x.Value).First().Key;
                ranked_mf_gbt.Add(max_mf_gbt);
                ranked_count_ByGBT.Add(Convert.ToDouble(max_mf_gbt.Length));
            }

            //plot conflict measure
            var plt_conflict = new ScottPlot.Plot(800, 400);
            plt_conflict.PlotScatter(time.ToArray(), conflictMeasure.ToArray());
            plt_conflict.SaveFig("conflict_measure.png");

            //plot hartley measure
            var plt_hartley = new ScottPlot.Plot(800, 400);
            plt_hartley.PlotScatter(time.ToArray(), hartleyMeasure.ToArray());
            plt_hartley.SaveFig("hartley_measure.png");

            //plot the decision rules based on the count of strings that were ranked highest in the evidences
            var plt_rank_belief = new ScottPlot.Plot(800, 400);
            plt_rank_belief.PlotScatter(time.ToArray(), ranked_count_ByBelief.ToArray());
            plt_rank_belief.SaveFig("rank_belief_count.png");

            var plt_rank_plausibility = new ScottPlot.Plot(800, 400);
            plt_rank_plausibility.PlotScatter(time.ToArray(), ranked_count_ByPlausibility.ToArray());
            plt_rank_plausibility.SaveFig("rank_plausibility_count.png");

            var plt_rank_pignistic = new ScottPlot.Plot(800, 400);
            plt_rank_pignistic.PlotScatter(time.ToArray(), ranked_count_ByPignistic.ToArray());
            plt_rank_pignistic.SaveFig("rank_pignistic_count.png");

            var plt_rank_gbt = new ScottPlot.Plot(800, 400);
            plt_rank_gbt.PlotScatter(time.ToArray(), ranked_count_ByGBT.ToArray());
            plt_rank_gbt.SaveFig("rank_gbt_count.png");

            // Now create a window for creating intrusions and causing attacks on the grid
            // Let say a cyber intrusion was performed on time [5,20] .... then [80,120]... then get control [200,250] 
            //create a list of tuple 
            List<Tuple<int, int>> physical_intrusions = new List<Tuple<int, int>>();
            Tuple<int, int> pevent1 = new Tuple<int, int>(20, 25); physical_intrusions.Add(pevent1);
            Tuple<int, int> pevent2 = new Tuple<int, int>(118, 125); physical_intrusions.Add(pevent2);
            Tuple<int, int> pevent3 = new Tuple<int, int>(230, 240); physical_intrusions.Add(pevent3);
            Tuple<int, int> pevent4 = new Tuple<int, int>(465, 490); physical_intrusions.Add(pevent4);
            Tuple<int, int> pevent5 = new Tuple<int, int>(860, 910); physical_intrusions.Add(pevent5);


            // generate physical data: voltage for a 3 bus system
            int nbus = 3;
            
            Random rand = new Random();
            List<List<double>> voltages = new List<List<double>>();
            for (int i = 0; i < nbus; i++)
            {
                List<double> volt_bus = new List<double>();
                for (int j = 0; j < sampleSize; j++)
                {
                    double start = 0.8;
                    double end = 1.2;
                    if (!IsInRange(j, physical_intrusions))
                    { end = 1.05; start = 0.95; }
                    var v = (rand.NextDouble() * Math.Abs(end - start)) + start;
                    volt_bus.Add(v);
                }
                voltages.Add(volt_bus);
            }

            //plot voltages using scottplot
            bool smoothen = true;
            if(smoothen)
            {
                voltages[0] = MovingAverage(10, voltages[0]);
                voltages[1] = MovingAverage(10, voltages[1]);
                voltages[2] = MovingAverage(10, voltages[2]);
            }


            var v1 = voltages[0].ToArray();
            var v2 = voltages[1].ToArray();
            var v3 = voltages[2].ToArray();
            var plt1 = new ScottPlot.Plot(800, 400);
            var plt2 = new ScottPlot.Plot(800, 400);
            var plt3 = new ScottPlot.Plot(800, 400);

            

            plt1.PlotScatter(time.ToArray(),v1);
            plt1.SaveFig("v1.png");

            plt2.PlotScatter(time.ToArray(), v2);
            plt2.SaveFig("v2.png");

            plt3.PlotScatter(time.ToArray(),v3);
            plt3.SaveFig("v3.png");

            // compute mass function for the physical values
            double t_low_lower_limit = 0.93;
            double t_high_lower_limit = 0.97;
            double t_low_higher_limit = 1.03;
            double t_high_higher_limit = 1.07;



            List<Dictionary<string, double>> phyList = new List<Dictionary<string, double>>();
            String[] phyKeys = new string[3] { "x","y","z" };
            for (int i = 0; i < sampleSize; i++)
            {
                Dictionary<string, double> mf_per_sample = new Dictionary<string, double>();
                for (int j = 0; j < nbus; j++)
                {
                    if(voltages[j][i] >= t_high_lower_limit && voltages[j][i] <= t_low_higher_limit)
                    {
                        mf_per_sample[phyKeys[j]] = 0.0;
                    }
                    else if(voltages[j][i] <= t_high_lower_limit && voltages[j][i] >= t_low_lower_limit)
                    {
                        mf_per_sample[phyKeys[j]] = 1.0 + Convert.ToDouble(Convert.ToDouble(t_low_lower_limit - voltages[j][i]) / Convert.ToDouble(t_high_lower_limit - t_low_lower_limit));
                    }
                    else if (voltages[j][i] <= t_high_higher_limit && voltages[j][i] >= t_low_higher_limit)
                    {
                        mf_per_sample[phyKeys[j]] = 1.0 + Convert.ToDouble(Convert.ToDouble(voltages[j][i] - t_high_higher_limit) / Convert.ToDouble(t_high_higher_limit - t_low_higher_limit));
                    }
                    else
                    {
                        mf_per_sample[phyKeys[j]] = 1.0;
                    }
                }
                if(IsAllBPANull(mf_per_sample))
                {
                    phyList.Add(new Dictionary<string, double> { { "", 1.0 } });
                }
                else
                {
                    phyList.Add(normalize(mf_per_sample));
                }                
            }

            List<string> cp_ranked_Hyp_ByBelief = new List<string>();
            List<double> cp_ranked_count_ByBelief = new List<double>();
            List<string> cp_ranked_Hyp_ByPlausibility = new List<string>();
            List<double> cp_ranked_count_ByPlausibility = new List<double>();

            List<string> cp_rankedPignistic = new List<string>();
            List<double> cp_ranked_count_ByPignistic = new List<double>();

            // construct mf using Generalized bayesian Theorem and store it here for every sample fused.
            List<Dictionary<string, double>> cp_gbt_sample_fused = new List<Dictionary<string, double>>();
            List<string> cp_ranked_mf_gbt = new List<string>();
            List<double> cp_ranked_count_ByGBT = new List<double>();

            List<double> cpConflictMeasure = new List<double>();
            List<double> cpHartleyMeasure = new List<double>();

            bool combinedCautiousCP = false;
            //test fusing cyber sensor with the physical sensor
            for (int i = 0; i < sampleSize; i++)
            {
                var s1_mf = sensor1_mf_list[i];
                var p1_mf = phyList[i];
                
                var fused = m.combinedDeterministically(s1_mf, p1_mf, true);
                if (combinedCautiousCP) fused = m.combineCautious(s1_mf, p1_mf,true);

                //conflict measure
                MassFunction t1 = new MassFunction(s1_mf);
                MassFunction t2 = new MassFunction(p1_mf);
                var conflict = m.conflict(t1, t2, 0);
                cpConflictMeasure.Add(conflict);

                //hartley measure
                var hartley = m.hartley_measure(fused);
                cpHartleyMeasure.Add(hartley);

                MassFunction newM = new MassFunction(fused);
                newM.ComputeFocalSet();
                newM.ComputeFrameOfDiscernment();
                var bel_func = newM.beliefFunction();
                var plaus_func = newM.plausFunction();

                //// test gbt 
                //var gbt_mf = newM.gbt(plaus_func, false, 0);
                //cp_gbt_sample_fused.Add(gbt_mf);

                // decision
                var max_hyp_bel = bel_func.OrderByDescending(x => x.Value).First().Key;
                cp_ranked_Hyp_ByBelief.Add(max_hyp_bel);
                cp_ranked_count_ByBelief.Add(Convert.ToDouble(max_hyp_bel.Length));

                var max_hyp_plaus = plaus_func.OrderByDescending(x => x.Value).First().Key;
                cp_ranked_Hyp_ByPlausibility.Add(max_hyp_plaus);
                cp_ranked_count_ByPlausibility.Add(Convert.ToDouble(max_hyp_plaus.Length));

                var pignisticFunction = newM.pignistic();
                if (pignisticFunction != null && pignisticFunction.Count != 0)
                {
                    var max_pignistic = pignisticFunction.OrderByDescending(x => x.Value).First().Key;
                    rankedPignistic.Add(max_pignistic);
                    cp_ranked_count_ByPignistic.Add(Convert.ToDouble(max_pignistic.Length));
                }
                else
                {
                    rankedPignistic.Add(max_hyp_plaus);
                    cp_ranked_count_ByPignistic.Add(Convert.ToDouble(max_hyp_plaus.Length));
                }

                //var max_mf_gbt = gbt_mf.OrderByDescending(x => x.Value).First().Key;
                //ranked_mf_gbt.Add(max_mf_gbt);
                //cp_ranked_count_ByGBT.Add(Convert.ToDouble(max_mf_gbt.Length));
            }

            //plot conflict measure
            var cp_plt_conflict = new ScottPlot.Plot(800, 400);
            cp_plt_conflict.PlotScatter(time.ToArray(), cpConflictMeasure.ToArray());
            cp_plt_conflict.SaveFig("cp_conflict_measure.png");

            //plot hartley measure
            var cp_plt_hartley = new ScottPlot.Plot(800, 400);
            cp_plt_hartley.PlotScatter(time.ToArray(), cpHartleyMeasure.ToArray());
            cp_plt_hartley.SaveFig("cp_hartley_measure.png");

            //plot the decision rules based on the count of strings that were ranked highest in the evidences
            var cp_plt_rank_belief = new ScottPlot.Plot(800, 400);
            cp_plt_rank_belief.PlotScatter(time.ToArray(), cp_ranked_count_ByBelief.ToArray());
            cp_plt_rank_belief.SaveFig("cp_rank_belief_count.png");

            var cp_plt_rank_plausibility = new ScottPlot.Plot(800, 400);
            cp_plt_rank_plausibility.PlotScatter(time.ToArray(), cp_ranked_count_ByPlausibility.ToArray());
            cp_plt_rank_plausibility.SaveFig("cp_rank_plausibility_count.png");

            var cp_plt_rank_pignistic = new ScottPlot.Plot(800, 400);
            cp_plt_rank_pignistic.PlotScatter(time.ToArray(), cp_ranked_count_ByPignistic.ToArray());
            cp_plt_rank_pignistic.SaveFig("cp_rank_pignistic_count.png");

            //var cp_plt_rank_gbt = new ScottPlot.Plot(800, 400);
            //cp_plt_rank_gbt.PlotScatter(time.ToArray(), cp_ranked_count_ByGBT.ToArray());
            //cp_plt_rank_gbt.SaveFig("cp_rank_gbt_count.png");



            Dictionary<string, double> testpf = new Dictionary<string, double> { { "b", 0.5 }, { "c", 0.8 } };
            var m_gbt = m.gbt(testpf, true, 20000);

            Dictionary<string, double> hyp = new Dictionary<string, double> { { "ab", 0.3 }, { "bc", 0.5 }, { "abc", 0.2 } };

            MassFunction m2 = new MassFunction(hyp);
            m2.ComputeFocalSet();
            m2.ComputeFrameOfDiscernment();
            m2.weight_function();

            var conf = m.conflict(m, m2, 1000);

        }

        public static List<Dictionary<string,double>> GenerateRandomMassFunction(int sampleSize)
        {
            List<Dictionary<string, double>> mf_list = new List<Dictionary<string, double>>();
            for (int i = 0; i < sampleSize; i++)
            {
                // generate random
                Dictionary<string, double> random_hyp = new Dictionary<string, double>();
                Random r = new Random();
                for (int j = 0; j < 4; j++)
                {
                    var randint = r.Next(1, 4);
                    var randomString = GenerateCoupon(randint);

                    random_hyp[randomString] = r.NextDouble();

                }
                var normalized = normalize(random_hyp);
                mf_list.Add(normalized);
            }
            return mf_list;
        }

        public static List<Dictionary<string, double>> GenerateRandomMassFunctionScenario(int sampleSize, List<Tuple<int,int>> attackTimes)
        {
            List<Dictionary<string, double>> mf_list = new List<Dictionary<string, double>>();
            for (int i = 0; i < sampleSize; i++)
            {
                // generate random
                Dictionary<string, double> random_hyp = new Dictionary<string, double>();
                Random r = new Random();
                for (int j = 0; j < 4; j++)
                {
                    var randint = r.Next(1, 4);
                    var randomString = GenerateCoupon(randint);

                    if(IsInRange(i,attackTimes)) random_hyp[randomString] = r.NextDouble();
                    else random_hyp[randomString] = 0.0;

                }
                if (IsAllBPANull(random_hyp))
                {
                    random_hyp.Clear();
                    random_hyp[""] = 1.0;
                }
                var normalized = normalize(random_hyp);
                mf_list.Add(normalized);
            }
            return mf_list;
        }

        public static bool IsInRange(int index, List<Tuple<int,int>> events)
        {
            foreach (var evt in events)
            {
                if (index > evt.Item1 && index < evt.Item2) return true;
            }
            return false;
        }

        public static bool IsAllBPANull(Dictionary<string,double> mf)
        {
            foreach (var value in mf.Values)
            {
                if (value != 0.0) return false;
            }
            return true;
        }

        public static string GenerateCoupon(int length)
        {
            Random random = new Random();
            string characters = "abcd";
            StringBuilder result = new StringBuilder(length);
            for (int i = 0; i < length; i++)
            {
                var item = characters[random.Next(characters.Length)];
                if(!result.ToString().Contains(item.ToString())) result.Append(item);
            }
            return result.ToString();
        }

        public static Dictionary<string, double> normalize(Dictionary<string, double> mf)
        {
            double mass_sum = mf.Sum(x => x.Value);
            var keyList = mf.Keys.ToList();
            if (mass_sum != 1.0)
            {
                foreach (var item in keyList)
                {
                    mf[item] = mf[item] / mass_sum;
                }
            }
            return mf;
        }

        public static List<double> MovingAverage(int period, List<double> Data)
        {
            var movingAverages = new double[Data.Count];
            var runningTotal = 0.0d;
            for (int i = 0; i < Data.Count; i++)
            {
                runningTotal += Data[i];
                if (i - period >= 0)
                {
                    var lost = Data[i - period];
                    runningTotal -= lost;
                    movingAverages[i] = runningTotal / period;
                }
                else movingAverages[i] = Data[i];
            }
            return movingAverages.ToList();
        }

    }
}
