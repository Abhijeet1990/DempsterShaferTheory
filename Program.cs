using System;
using System.Collections.Generic;
using static SwagShopLibrary.DSTheory;


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


            Dictionary<string, double> hyp_val2 = new Dictionary<string, double> { { "ab", 0.4 }, { "bc", 0.4 }, { "a", 0.2 }, { "ad", 0.0 } };
            var combined = m.combinedDeterministically(hyp_val2, hyp_val, false); // conjunctive combination
            var combinedSample = m.combinedDirectSampling(hyp_val2, hyp_val, 1000, false); // conjunctive combination
            var combinedImportanceSample = m.combinedImportanceSampling(hyp_val2, hyp_val, 1000);
            var samples = m.sample(1000, true, mf_pf);


            Dictionary<string, double> testpf = new Dictionary<string, double> { { "b", 0.5 }, { "c", 0.8 } };
            var m_gbt = m.gbt(testpf, true, 20000);

            Dictionary<string, double> hyp = new Dictionary<string, double> { { "ab", 0.3 }, { "bc", 0.5 }, { "abc", 0.2 } };

            MassFunction m2 = new MassFunction(hyp);
            m2.ComputeFocalSet();
            m2.ComputeFrameOfDiscernment();
            m2.weight_function();

            Dictionary<string, double> mcc1 = new Dictionary<string, double> { { "ab", 0.3 }, { "bc", 0.5 }, { "abc", 0.2 } };
            Dictionary<string, double> mcc2 = new Dictionary<string, double> { { "b",0.3}, {"bc",0.4}, {"abc",0.3} };
            m.combineCautious(mcc1, mcc2);
        }
    }
}
