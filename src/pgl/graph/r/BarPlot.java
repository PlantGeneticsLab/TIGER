/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package pgl.graph.r;

import java.io.File;
import rcaller.RCaller;
import rcaller.RCode;

/**
 * Plot bar plot
 * @author Fei Lu
 */
public class BarPlot extends RgraphicsAbstract {
    double[] value = null;
    String[] name = null;
    
    public BarPlot (double[] value, String[] name) {
        this.value = value;
        this.name = name;
    }
    
    public BarPlot (int[] intValue, String[] name) {
        value = new double[intValue.length];
        for (int i = 0; i < value.length; i++) value[i] = intValue[i];
        this.name = name;
    }
    
    @Override
    public void showGraph() {
        RCaller caller = new RCaller();
        caller.setRscriptExecutable(this.rPath);
        RCode rCode = new RCode();
        rCode.addDoubleArray("x", value);
        rCode.addStringArray("name", name);
        try {
            File plotFile = rCode.startPlot();
            rCode.addRCode(margin);
            rCode.addRCode(this.getPlotStatement());
            caller.setRCode(rCode);
            caller.runOnly();
            rCode.endPlot();
            rCode.showPlot(plotFile);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        String s = rCode.getCode().toString();
        System.out.println(s);
    }
    
    @Override
    public void saveGraph(String outfileS) {
        outfileS = outfileS.replaceAll("\\\\", "/");
        RCaller caller = new RCaller();
        caller.setRscriptExecutable(this.rPath);
        RCode rCode = new RCode();
        rCode.addDoubleArray("x", value);
        rCode.addStringArray("name", name);
        rCode.addRCode("pdf(\""+outfileS+"\")");
        rCode.addRCode(margin);
        rCode.addRCode(this.getPlotStatement());
        rCode.addRCode("dev.off()");
        caller.setRCode(rCode);
        caller.runOnly();
    }

    @Override
    public String getPlotStatement() {
        StringBuilder sb = new StringBuilder ();
        sb.append("barplot(x, names.arg = name, ").append("col=\"").append(defaultColor).append("\", ");
        sb.append("main=\"").append(title).append("\", ");
        sb.append("ylab=\"").append(yLab).append("\", ");
        sb.append("cex.main=").append(titleSize).append(", ");
        sb.append("cex.lab=").append(labSize).append(", ");
        sb.append("cex.axis=").append(axisSize).append(", ");
        sb.append("font.main=").append(titleFont).append(",");
        if (sb.charAt(sb.length()-1) == ',') {
            sb.deleteCharAt(sb.length()-1);
        }
        sb.append(")");
        String statement = sb.toString();
        //System.out.println(statement);
        return statement;
    }
    
}
