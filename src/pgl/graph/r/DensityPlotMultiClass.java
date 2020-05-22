/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package pgl.graph.r;

import java.io.File;
import java.util.ArrayList;
import rcaller.RCaller;
import rcaller.RCode;

/**
 * Plot density plot of multiple classes in one graph, using ggplot2 package
 * @author Fei Lu
 */
public class DensityPlotMultiClass extends Rgraphics {
    ArrayList<double[]> valueList = null;
    String[] variableNames = null;
    
    public DensityPlotMultiClass (String[] variableNames) {
        this.variableNames = variableNames;
    }
    
    public DensityPlotMultiClass (double[][] values, String[] variableNames) {
        valueList = new ArrayList();
        for (int i = 0; i < values.length; i++) {
            valueList.add(values[i]);
        }
        this.variableNames = variableNames;
    }
    
    public DensityPlotMultiClass (int[][] values, String[] variableNames) {
        valueList = new ArrayList();
        for (int i = 0; i < values.length; i++) {
            double[] temp = new double[values[i].length];
            for (int j = 0; j < temp.length; j++) {
                temp[j] = values[i][j];
            }
            valueList.add(temp);
        }
        this.variableNames = variableNames;
    }
    
    public DensityPlotMultiClass (ArrayList<double[]> valueList, String[] variableNames) {
        this.valueList = valueList;
        this.variableNames = variableNames;
    }
    
    public void setIntValues (ArrayList<int[]> intValueList) {
        valueList = new ArrayList();
        double[] a;
        for (int i = 0; i < intValueList.size(); i++) {
            a = new double[intValueList.get(i).length];
            for (int j = 0; j < a.length; j++) {
                a[j] = intValueList.get(i)[j];
            }
            valueList.add(a);
        }
    }
    
    public void setDoubleValues (ArrayList<double[]> doubleValueList) {
        this.valueList = doubleValueList;
    }

    private String[] getFactorNames () {
        int sum = 0;
        for (int i = 0; i < this.valueList.size(); i++) {
            sum+=valueList.get(i).length;
        }
        String[] factorNames = new String[sum];
        int cnt = 0;
        for (int i = 0; i < this.valueList.size(); i++) {
            for (int j = 0; j < valueList.get(i).length; j++) {
                factorNames[cnt] = this.variableNames[i];
                cnt++;
            }
        }
        return factorNames;
    }
    
    @Override
    public void showGraph() {
        RCaller caller = new RCaller();
        caller.setRscriptExecutable(this.rPath);
        RCode rCode = new RCode();
        rCode.addRCode("library(ggplot2)");
        for (int i = 0; i < this.valueList.size(); i++) {
            rCode.addDoubleArray("x"+String.valueOf(i), this.valueList.get(i));
        }
        rCode.addStringArray("name", this.getFactorNames());
        String dfStr = "df <- data.frame(Class = factor(name), x = c(";
        for (int i = 0; i < this.valueList.size(); i++) {
            dfStr = dfStr + "x"+String.valueOf(i)+",";
        }
        dfStr = dfStr.substring(0, dfStr.length()-1) + "))";
        rCode.addRCode(dfStr);
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
    }

    @Override
    public void saveGraph(String outfileS) {
        outfileS = outfileS.replaceAll("\\\\", "/");
        RCaller caller = new RCaller();
        caller.setRscriptExecutable(this.rPath);
        RCode rCode = new RCode();
        rCode.addRCode("library(ggplot2)");
        for (int i = 0; i < this.valueList.size(); i++) {
            rCode.addDoubleArray("x"+String.valueOf(i), this.valueList.get(i));
        }
        rCode.addStringArray("name", this.getFactorNames());
        String dfStr = "df <- data.frame(Class = factor(name), x = c(";
        for (int i = 0; i < this.valueList.size(); i++) {
            dfStr = dfStr + "x"+String.valueOf(i)+",";
        }
        dfStr = dfStr.substring(0, dfStr.length()-1) + "))";
        rCode.addRCode(dfStr);
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
        sb.append("qplot(x, data = df, geom = \"density\", fill = Class, alpha = I(.5),");
        sb.append("main=\"").append(title).append("\",");
        sb.append("xlab=\"").append(xLab).append("\",");
        sb.append("ylab=\"").append(yLab).append("\",");
        sb.append("cex.main=").append(titleSize).append(",");
        sb.append("cex.lab=").append(labSize).append(",");
        sb.append("cex.axis=").append(axisSize).append(",");
        sb.append("font.main=").append(titleFont).append(",");
        if (!Double.isNaN(this.xLimLow)) {
            sb.append("xlim = c(").append(this.xLimLow).append(",").append(this.xLimHigh).append("),");
        }
        if (!Double.isNaN(this.yLimLow)) {
            sb.append("ylim = c(").append(this.yLimLow).append(",").append(this.yLimHigh).append("),");
        }
        if (sb.charAt(sb.length()-1) == ',') {
            sb.deleteCharAt(sb.length()-1);
        }
        sb.append(")");
        String statement = sb.toString();
        //System.out.println(statement);
        return statement;
    }
}
