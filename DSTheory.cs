using Combinatorics.Collections;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SwagShopLibrary
{
    public class DSTheory
    {
        public class MassFunction
        {
            //A Dempster-Shafer mass function(basic probability assignment) based on a dictionary.
            //Both normalized and unnormalized mass functions are supported.
            //The underlying frame of discernment is assumed to be discrete.
            //Hypotheses and their associated mass values can be added/changed/removed using the standard dictionary methods.
            //Each hypothesis can be an arbitrary sequence which is automatically converted to a 'frozenset', meaning its elements must be hashable.
            public Dictionary<string, double> hyp_value;
            public Dictionary<string, double> focal_set = new Dictionary<string, double>();
            public string fullList = "";
            public List<string> frameOfDiscernment = new List<string>();

            public MassFunction(Dictionary<string, double> _hyp_value)
            {
                hyp_value = _hyp_value;
            }
           
            public void ComputeFocalSet()
            {
                foreach (var item in hyp_value.Keys)
                {
                    if(hyp_value[item]!=0)
                    {
                        if(!focal_set.ContainsKey(item))
                        {
                            focal_set.Add(item, hyp_value[item]);
                        }
                    }
                }
            }

            public void ComputeFrameOfDiscernment()
            {
                string toSend = "";
                foreach (var item in focal_set.Keys)
                {
                    foreach(var c in item)
                    {
                        if(!toSend.Contains(c))
                        {
                            toSend += c;
                        }
                    }
                }
                if (fullList == "") fullList = toSend;
                var list = CalculateCombinationsList(toSend);
                foreach(var i in list)
                {
                    if(!frameOfDiscernment.Contains(i))
                    {
                        frameOfDiscernment.Add(i);
                    }
                }
                
            }

            public List<string> CalculateCombinationsList(string str)
            {
                var results = Enumerable.Range(0, 1 << str.Length)
                .Select(e => string.Join(string.Empty, Enumerable.Range(0, str.Length).Select(b => (e & (1 << b)) == 0 ? (char?)null : str[b])));
                return results.ToList();
            }

            public double belief(string hyp)
            {
                double belief = 0.0;
                foreach (var i in hyp_value.Keys)
                {
                    if (hyp.Contains(i))  // y(key in mass function) subset of x(hyp)
                    {
                        belief += hyp_value[i];
                    }

                }
                return belief;
            }
            public Dictionary<string, double> beliefFunction()
            {
                Dictionary<string, double> beliefFunction = new Dictionary<string, double>();
                foreach (var i in frameOfDiscernment)
                {
                    if(!beliefFunction.ContainsKey(i))
                    {
                        beliefFunction.Add(i, belief(i));
                    }
                }
                return beliefFunction;
            }

            public double plausibility(string hyp)
            {
                double plausibility = 0.0;
                foreach(var i in hyp_value.Keys)
                {
                    foreach(char c in hyp)
                    {
                        if (i.Contains(c))
                        {
                            plausibility += hyp_value[i];
                            break;
                        }
                    }
                }
                return plausibility;
            }
            public Dictionary<string, double> plausFunction()
            {
                Dictionary<string, double> plausFunction = new Dictionary<string, double>();
                foreach (var i in frameOfDiscernment)
                {
                    if (!plausFunction.ContainsKey(i))
                    {
                        plausFunction.Add(i, plausibility(i));
                    }
                }
                return plausFunction;
            }

            public double commonality(string hyp)
            {
                double commonality = 0.0;
                foreach (var i in hyp_value.Keys)
                {
                    if (i.Contains(hyp))  // x(hyp) subset of y(key in mass function)
                    {
                        commonality += hyp_value[i];
                    }

                }
                return commonality;
            }
            public Dictionary<string, double> commFunction()
            {
                Dictionary<string, double> commFunction = new Dictionary<string, double>();
                foreach (var i in frameOfDiscernment)
                {
                    if (!commFunction.ContainsKey(i))
                    {
                        commFunction.Add(i, commonality(i));
                    }
                }
                return commFunction;
            }

            // compute mass function from belief function
            public Dictionary<string, double> fromBelief(Dictionary<string,double> belFunc )
            {
                Dictionary<string, double> massFromBelief = new Dictionary<string, double>();
                //foreach (var item in frameOfDiscernment)
                foreach (var item in belFunc.Keys)
                    {
                    int aCount = item.Length;
                    double sum = 0.0;
                    foreach (var key in CalculateCombinationsList(item))
                    {
                        int bCount = key.Length;
                        //if (!item.Contains(key))
                        //{
                            sum += Math.Pow(-1, aCount - bCount) * belFunc[key];
                        //}
                    }
                    if(!massFromBelief.ContainsKey(item) && sum>0)
                    {
                        massFromBelief[item] = sum;
                    }
                }
                return massFromBelief;
            }

            // compute mass function from plausibility function
            public Dictionary<string, double> fromPlausibility(Dictionary<string, double> plaFunc)
            {
                Dictionary<string, double> belief = new Dictionary<string, double>();
                
                // get the key of the maximum value in the plausible function
                var pl_max = plaFunc.Aggregate((l, r) => l.Value > r.Value ? l : r).Key;
                var bel_theta = plaFunc[pl_max];
                foreach(var item in plaFunc.Keys)
                {
                    var dummy = fullList;
                    foreach (var j in item)
                    {
                        int idx = dummy.IndexOf(j);
                        if (idx != -1) dummy=dummy.Remove(idx, 1);
                    }
                    if (!belief.ContainsKey(dummy))
                    {
                        belief[dummy] = bel_theta - plaFunc[item];
                    }
                }
                return fromBelief(belief);
               
            }

            //compute mass function from commonality function

            public string Filter(string str, List<char> charsToRemove)
            {
                foreach (char c in charsToRemove) str = str.Replace(c.ToString(), String.Empty);
                return str;
            }
            public Dictionary<string,double> fromCommonality(Dictionary<string,double> commFunc)
            {
                Dictionary<string, double> massFromComm = new Dictionary<string, double>();

                // get the key of the maximum value in the plausible function
                //commFunc.Remove("");
                var q_max = commFunc.Aggregate((l, r) => l.Value > r.Value ? l : r).Key;
                foreach (var h1 in commFunc.Keys)
                {
                    //int aCount = h1.Length;
                    double sum = 0.0;
                    List<char> h1_list = new List<char>(h1);
                    var remaining = Filter(q_max, h1_list);
                    foreach (var h2 in CalculateCombinationsList(remaining))
                    {
                        //int bCount = h2.Length;
                        List<char> h2_list = new List<char>(h2);
                        var union = h1_list.Union(h2_list).ToList();

                        foreach (var key in commFunc.Keys)
                        {
                            var chars = new List<char>(key);
                            var check = union.All(chars.Contains) && union.Count == chars.Count;
                            if (check) sum += Math.Pow(-1, (Filter(h2,h1_list).Length)) * commFunc[key]; 
                        }
                    }
                    if (!massFromComm.ContainsKey(h1) && sum > 0) massFromComm[h1] = sum;

                }

                return massFromComm;
            }

            public Dictionary<string,double> combine(Dictionary<string, double> mf1, Dictionary<string, double> mf2, bool normalization, int sample_count, bool importance_sampling, bool isDisjunctive)
            {
                Dictionary<string, double> combined = new Dictionary<string, double>();
                if(sample_count == 0)
                {
                    combined = combinedDeterministically(mf1, mf2, isDisjunctive);
                }
                else if(importance_sampling)
                {
                    combined = combinedImportanceSampling(mf1, mf2, sample_count);
                }
                else
                {
                    combined = combinedDirectSampling(mf1, mf2, sample_count, isDisjunctive);
                }

                if (normalization)
                {
                    return normalize(combined);
                }
                else return combined;

            }

            // Deterministic combination of mass function
            public Dictionary<string, double> combinedDeterministically(Dictionary<string,double> mf1, Dictionary<string,double> mf2, bool isDisjunctive)
            {
                Dictionary<string, double> combined = new Dictionary<string, double>();
                foreach (var item in mf1.Keys)
                {
                    foreach (var item2 in mf2.Keys)
                    {
                        var mixed = item.Intersect(item2);
                        var m = String.Concat(mixed);
                        if(isDisjunctive)
                        {
                            mixed = item.Union(item2);
                            m = new string(String.Concat(mixed).ToCharArray().Distinct().ToArray());
                        }


                        //if abc was there but say bac came then they must add. rather than create a new key
                        foreach (var i in combined.Keys) if (checkEquality(m,i)) m = i;
                        
                        if (!combined.ContainsKey(m))
                        {
                            combined[m]= mf1[item] * mf2[item2];
                        }
                        else
                        {
                            combined[m]+= mf1[item] * mf2[item2];
                        }
                            
                    }
                }

                // normalize to remove the ignorance
                combined = normalizeRemoveIgnorance(combined);

                return combined;
            }

            public Dictionary<string, double> combinedDirectSampling(Dictionary<string, double> mf1, Dictionary<string, double> mf2, int sample_count,bool isDisjunctive)
            {
                Dictionary<string, double> combined = new Dictionary<string, double>();
                var samples1 = sample(sample_count, true, mf1);
                var samples2 = sample(sample_count, true, mf2);
                for (int i = 0; i < sample_count; i++)
                {
                    var mixed = samples1[i].Intersect(samples2[i]);
                    var m = String.Concat(mixed);
                    if (isDisjunctive)
                    {
                        mixed = samples1[i].Union(samples2[i]);
                        m = new string(String.Concat(mixed).ToCharArray().Distinct().ToArray());
                    }

                    if (!combined.ContainsKey(m))
                    {
                        combined[m] = (double)1/(double)sample_count;
                    }
                    else
                    {
                        combined[m] += (double)1 / (double)sample_count;
                    }

                }
                return combined;
            }

            public Dictionary<string, double> combinedImportanceSampling(Dictionary<string, double> mf1, Dictionary<string, double> mf2, int sample_count)
            {
                Dictionary<string, double> combined = new Dictionary<string, double>();
                var sampled = sampleAsDict(sample_count, true, mf1);
                foreach (var item in sampled.Keys)
                {
                    var weight = plausibility(item);
                    Dictionary<string, double> fix = new Dictionary<string, double> { { item, 1.0 } };
                    foreach (var item2 in sample(sampled[item],true,fix))
                    {
                        if (!combined.ContainsKey(item2))
                        {
                            combined[item2] = weight;
                        }
                        else combined[item2] += weight;
                    }
                }
                return combined;
            }

            public List<string> sample(int n,bool quantization, Dictionary<string,double> mf)
            {
                List<string> samples = new List<string>();
                double mass_sum = mf.Sum(x => x.Value);

                if(quantization)
                {
                    Dictionary<string, double> remainders = new Dictionary<string, double>();
                    var remaining_sample_count = n;
                    foreach (var item in mf.Keys)
                    {
                        double fraction = n * mf[item] / mass_sum;
                        int quotient = Convert.ToInt32(fraction);
                        for (int i = 0; i < quotient; i++)
                        {
                            samples.Add(item);
                        }
                        remainders.Add(item, fraction - quotient);
                        remaining_sample_count -= quotient;
                    }
                    // sort by descending order the dictionary
                    var x = from pair in remainders
                                 orderby pair.Value descending
                                 select pair;
                    remainders = x.ToDictionary(pair => pair.Key, pair => pair.Value);

                    foreach (var item in remainders.Keys.Take(remaining_sample_count))
                    {
                        samples.Add(item);
                    }
                }
                return samples;
            }

            public Dictionary<string,int> sampleAsDict(int n, bool quantization, Dictionary<string, double> mf)
            {
                Dictionary<string,int> samples =new Dictionary<string, int>();
                double mass_sum = mf.Sum(x => x.Value);

                if (quantization)
                {
                    Dictionary<string, double> remainders = new Dictionary<string, double>();
                    var remaining_sample_count = n;
                    foreach (var item in mf.Keys)
                    {
                        double fraction = n * mf[item] / mass_sum;
                        int quotient = Convert.ToInt32(fraction);
                        samples[item] = quotient;                       
                        remainders.Add(item, fraction - quotient);
                        remaining_sample_count -= quotient;
                    }
                    // sort by descending order the dictionary
                    var x = from pair in remainders
                            orderby pair.Value descending
                            select pair;
                    remainders = x.ToDictionary(pair => pair.Key, pair => pair.Value);

                    foreach (var item in remainders.Keys.Take(remaining_sample_count))
                    {
                        samples[item]+=1;
                    }
                }
                return samples;
            }

            public Dictionary<string,double> normalize(Dictionary<string,double> mf)
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

            public Dictionary<string, double> normalizeRemoveIgnorance(Dictionary<string, double> mf)
            {
                double mass_sum = mf.Where(x=>x.Key!="").Sum(x => x.Value);
                var keyList = mf.Keys.ToList();
                if (mass_sum != 1.0)
                {
                    foreach (var item in keyList)
                    {
                        mf[item] = mf[item] / mass_sum;
                    }
                }
                if (mf.ContainsKey("")) mf.Remove("");
                return mf;
            }
            // Got to fix this one

            /*Constructs a mass function using the generalized Bayesian theorem.
        For more information, see Smets. 1993. Belief functions: 
        Paper: The disjunctive rule of combination and the generalized Bayesian theorem. International Journal of Approximate Reasoning. 
        
        'likelihoods' specifies the conditional plausibilities for a set of singleton hypotheses.
        It can either be a dictionary mapping singleton hypotheses to plausibilities or an iterable
        containing tuples consisting of a singleton hypothesis and a corresponding plausibility value.     
        'normalization' determines whether the resulting mass function is normalized, i.e., whether m({}) == 0.       
        If 'sample_count' is not None, the true mass function is approximated using the specified number of samples.*/
            public Dictionary<string, double> gbt(Dictionary<string, double>  likelihoods,bool normalization, int sample_count)
            {
                Dictionary<string, double> m = new Dictionary<string, double>();
                // dict.Where(item => item.Value.position.Equals(newPosition)).Select(item => item.Key)
                List<string> ones = new List<string>();
                foreach (var item in likelihoods.Keys)
                {
                    var k = likelihoods[item];
                    if (likelihoods[item]>=1)
                    {
                        ones.Add(item);
                    }
                }
                foreach (var item in likelihoods.Where(x=>x.Value>=1).ToList())
                {
                    likelihoods.Remove(item.Key);
                }
                //likelihoods = likelihoods.Where(x => x.Value <= 1).ToDictionary(pair => pair.Key, pair => pair.Value);
                var one = new string(String.Concat(ones).ToCharArray().Distinct().ToArray());
                if (sample_count == 0)
                {
                    traverse(m, likelihoods, one, 0, "", 1.0);
                    if (normalization)
                    {
                        normalize(m);
                    }
                }
                else // Monte Carlo Simulation
                {
                    double empty_mass = 1.0;
                    if(normalization)
                    {
                        foreach (var item in likelihoods.Keys)
                        {
                            empty_mass *= likelihoods[item];
                        }
                    }
                    for (int i = 0; i < sample_count; i++)
                    {
                        List<double> rv = new List<double>();
                        Random rnd = new Random();
                        foreach (var item in likelihoods.Keys)
                        {
                            rv.Add(rnd.NextDouble());
                        }
                        double subtree_mass = 1.0;
                        List<string> hyp = new List<string>();
                        hyp.AddRange(ones);
                        foreach (var item in likelihoods.Keys)
                        {
                            var l = likelihoods[item];
                            var p_t = l * subtree_mass;
                            var p_f = (1.0 - l) * subtree_mass;
                            if(normalization && hyp.Count==0)
                            {
                                p_f -= empty_mass;
                            }
                            int ind = Array.IndexOf(likelihoods.Keys.ToArray(), item);
                            if (p_t > rv[ind]*(p_t+p_f))
                            {
                                hyp.Add(item);
                            }
                            else
                            {
                                subtree_mass *= (1 - l);
                            }

                        }
                        var index_mf = new string(String.Concat(hyp).ToCharArray().Distinct().ToArray());
                        if (m.ContainsKey(index_mf))
                        {
                            m[index_mf] += (double)1 / (double)sample_count;
                        }
                        else
                        {
                            m[index_mf] = (double)1 / (double)sample_count;
                        }
                    }

                }
                
                return m;

            }

            private void traverse(Dictionary<string,double> m, Dictionary<string, double> likelihoods,string ones, int index, string hyp, double mass)
            {
                if (index == likelihoods.Count)
                {
                    var unique = new string(String.Concat(hyp+ones).ToCharArray().Distinct().ToArray());
                    var found = false;
                    foreach (var key in m.Keys)
                    {
                        var keychars = new List<char>(key);
                        var check = checkCharList(new List<char>(unique), keychars);
                        if (check)
                        {
                            m[key] += mass;
                            found = true;
                            break;
                        } 
                    }
                    if(!found) m[unique] = mass;
                }
                else
                {
                    traverse(m, likelihoods, ones, index + 1, hyp + likelihoods.Keys.ElementAt(index), mass * likelihoods.Values.ElementAt(index));
                    traverse(m, likelihoods, ones, index + 1, hyp, mass * (1 - likelihoods.Values.ElementAt(index)));
                }
            }

            public bool checkCharList(List<char> a, List<char> b)
            {
                return a.All(b.Contains) && a.Count == b.Count;
            }

            /*Combines the mass function with another mass function using the cautious rule and returns the combination as a new mass function.
        
        For more details, see:
        T. Denoeux (2008), "Conjunctive and disjunctive combination of belief functions induced by nondistinct bodies of evidence",
        Artificial Intelligence 172, 234-264.*/
            public Dictionary<string,double> combineCautious(Dictionary<string, double> mf1, Dictionary<string, double> mf2, bool disjunctive)
            {
                MassFunction m1 = new MassFunction(mf1);
                m1.ComputeFocalSet();
                m1.ComputeFrameOfDiscernment();
                
                MassFunction m2 = new MassFunction(mf2);
                m2.ComputeFocalSet();
                m2.ComputeFrameOfDiscernment();
                

                string all = "", theta = "";
                foreach (var item in mf1.Keys)
                {
                   all += item;
                   theta = new string(String.Concat(all).ToCharArray().Distinct().ToArray());
                }
                var w1 = m1.weight_function();
                var w2 = m2.weight_function();
                Dictionary<string, double> wmin = new Dictionary<string, double>();
                foreach (var item in w1.Keys)
                {
                    var replacedItem = item;
                    // if the key is not present...then change the ordering of the key based on the key of the first mass function
                    if (!w2.ContainsKey(item))
                    {
                        List<string> d = new List<string>();
                        permute(item.ToCharArray(), 0, item.Length - 1,ref d);
                        foreach (var i in w2.Keys)
                        {
                            if (d.Contains(i)) replacedItem = i;
                        }
                    }
                    if(w2.ContainsKey(replacedItem))
                    {
                        var x = w1[item] < w2[replacedItem] ? w1[item] : w2[replacedItem];
                        wmin.Add(item, x);
                    }
                    
                }
                Dictionary<string, double> one = new Dictionary<string, double> { { theta,1.0} };
                MassFunction m = new MassFunction(one);
                m.ComputeFocalSet();
                m.ComputeFrameOfDiscernment();
                

                foreach (var item in wmin.Keys)
                {
                    MassFunction m_simple = new MassFunction( new Dictionary<string, double> { { theta, wmin[item] }, { item, (1.0 - wmin[item]) } });
                    m_simple.ComputeFocalSet();
                    m_simple.ComputeFrameOfDiscernment();
                    if (disjunctive)
                    { 
                        m.hyp_value = combine(m.hyp_value, m_simple.hyp_value, false, 0, false, true); 
                    }
                    else
                    {
                        m.hyp_value = combine(m.hyp_value, m_simple.hyp_value, false, 0, false, false);
                    }

                }

                return m.hyp_value;
            }


            //public Dictionary<string,double> markov()
            //Computes the weight function corresponding to this mass function. 
            // The weights are computed based on the commonality function. 
            public Dictionary<string,double> weight_function()
            {
                Dictionary<string, double> weights = new Dictionary<string, double>();
                var q = commFunction();
                List<string> theta = new List<string>();
                
                for (int i = 0; i < frameOfDiscernment.Count; i++)
                {
                    theta.Add(frameOfDiscernment[i]);
                    
                }
                int maxLen = theta.Aggregate("", (max, cur) => max.Length > cur.Length ? max : cur).Length;
                for (int i = 0; i <theta.Count; i++)
                {
                    List<string> temp = new List<string>();
                    for (int j = 0; j < theta.Count; j++)
                    {
                        temp.Add(theta[j]);
                    }
                    if(theta[i].Length < maxLen)
                    {
                        //likelihoods.Where(x=>x.Value>=1)
                        if(theta[i]!="") temp.Remove(theta[i]);
                        List<string> sets = new List<string>();
                        foreach (var item in temp)
                        {
                            var h = new string(String.Concat(item + theta[i]).ToCharArray().Distinct().ToArray());
                            if (q.ContainsKey(h) && !sets.Contains(h))
                                sets.Add(h);

                        }
                        double q_even = 1.0;
                        double q_odd = 1.0;
                        for (int k = 0; k < sets.Count; k++)
                        {
                            if (sets[k].Length % 2 == 0) q_even *= q[sets[k]];
                            else q_odd *= q[sets[k]];
                        }
                        if(theta[i].Length%2==0)
                        {                
                            weights[theta[i]] = q_odd / q_even;
                        }
                        else
                        {
                            weights[theta[i]] = q_even/ q_odd;                       
                        }
                        if (q_odd == 0 || q_even == 0) weights[theta[i]] = 1.0;
                    }
                }
                return weights;
            }

            static void permute(char[] arry, int i, int n, ref List<string> d)
            {
                int j;
                if (i == n)
                {
                    string f = new string(arry.Distinct().ToArray());
                    d.Add(f);
                }
                else
                {
                    for (j = i; j <= n; j++)
                    {
                        swap(ref arry[i], ref arry[j]);
                        permute(arry, i + 1, n, ref d);
                        swap(ref arry[i], ref arry[j]); //backtrack
                    }
                }
            }

            static bool checkEquality(string a, string b)
            {
                if (a.Length == b.Length)
                {
                    foreach (char c in a)
                    {
                        if (b.Contains(c)) continue;
                        else return false;
                    }
                }
                else return false;
                return true;
            }

            static void swap(ref char a, ref char b)
            {
                char tmp;
                tmp = a;
                a = b;
                b = tmp;
            }

            //Calculates the weight of conflict between two or more mass functions.
            public double conflict(MassFunction m1, MassFunction m2, int sample_count)
            {
                var combined = combine(m1.hyp_value, m2.hyp_value, false,sample_count, false, false);
                var sum = combined.Sum(x => x.Value);
                if (sum == 0.0) return 10000000.0; // return some infinite number
                else return -Math.Log(sum);              
            }

            // Computes the local conflict measure.
            // returns 0 for unnormalized mass function
            public double local_conflict()
            {
                double sum = 0.0;
                foreach (var key in this.hyp_value.Keys)
                {
                    var value = this.hyp_value[key];
                    if (key == "" && this.hyp_value[key] > 0.0) return 0.0;                   
                    else sum += value * Math.Log(key.Length / value, 2);
                }
                return sum;
            }

            //  Removes all non-focal (0 mass) hypotheses in-place.
            public void prune()
            {
                var item = this.hyp_value.First(kvp => kvp.Value == 0.0);

                this.hyp_value.Remove(item.Key);
            }

            // computes the p-norm between two mass function
            public double norm(MassFunction m1, MassFunction m2, int p)
            {
                double sum = 0.0;
                foreach (var key in m1.hyp_value.Keys)
                {
                    if(m2.hyp_value.ContainsKey(key)) sum += Math.Pow(m1.hyp_value[key] - m2.hyp_value[key], p);
                }
                foreach (var key in m2.hyp_value.Keys)
                {
                    if (!m1.hyp_value.ContainsKey(key))
                    {
                        sum += Math.Pow(m2.hyp_value[key], p);
                    }
                }

                return Math.Pow(sum,1.0/p);
            }
            //Computes the Hartley-like measure in order to quantify the amount of imprecision.

            //For more information, see:
        //G.J.Klir(1999), "Uncertainty and information measures for imprecise probabilities: An overview",
        //International Symposium on Imprecise Probabilities and Their Applications.
            public double hartley_measure(Dictionary<string,double> hyp_value )
            {
                double sum = 0.0;
                foreach (var key in hyp_value.Keys)
                {
                    if (key.Length == 0) sum += 0;
                    else sum += hyp_value[key] * Math.Log(key.Length, 2.0);
                }
                return sum;
            }

            //Computes the pignistic transformation and returns it as a new mass function consisting only of singletons.
            public Dictionary<string,double> pignistic()
            {
                Dictionary<string, double> p = new Dictionary<string, double>();
                foreach (var key in this.hyp_value.Keys)
                {
                    if (this.hyp_value[key] > 0.0)
                    {
                        var charList = new List<char>(key);
                        foreach ( var character in charList)
                        {
                            if (p.ContainsKey(character.ToString()))
                            {
                                p[character.ToString()] += this.hyp_value[key] / charList.Count;
                            }
                            else
                            {
                                p.Add(character.ToString(), this.hyp_value[key] / charList.Count);
                            }
                        }
                    }

                }

                return normalize(p);
            }
        }
    }
}
