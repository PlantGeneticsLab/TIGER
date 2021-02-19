package pgl.app.fastCall2;

import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import pgl.AppUtils;
import pgl.PGLConstraints;
import pgl.infra.dna.FastaBit;
import pgl.infra.utils.*;
import java.io.BufferedReader;
import java.io.DataOutputStream;
import java.io.File;
import java.io.InputStreamReader;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.LongAdder;

public class FastCall2 {
    //Current step ID of the pipeline
    int step = Integer.MIN_VALUE;
    //genome block for individual ing
    static int binSize = 50000;
    //+, -, A, C, G, T,
    static final byte[] pileupAlleleAscIIs = {43, 45, 65, 67, 71, 84};

    public FastCall2 (String parameterFileS) {
        this.runSteps(parameterFileS);
    }

    private void runSteps(String parameterFileS) {
        Dyad<List<String>, List<String>> d = AppUtils.getParameterList(parameterFileS);
        List<String> pLineList = d.getFirstElement();
        List<String> sLineList = d.getSecondElement();
        this.step = Integer.parseInt(sLineList.get(0).split("\\s+")[1]);
        if (step == 1) {
            new DiscoverVariation(pLineList);
            //this.parseParametersStep1(pLineList);
        }
        else if (step == 2) {

        }
    }

    static byte getCodedDepth (int depth) {
        int dep = -128+depth;
        return (byte)dep;
    }

    static int getDecodeDepth (byte dep) {
        return dep+128;
    }
}
