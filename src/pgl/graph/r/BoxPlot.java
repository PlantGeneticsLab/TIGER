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
 * Plot box plot
 * @author Fei Lu
 */
public class BoxPlot extends RgraphicsAbstract {
    ArrayList<double[]> valueList = null;
    String[] variableNames = null;
    String[] colors = null;
    boolean ifXVertical = false;
    public BoxPlot (String[] variableNames) {
        this.variableNames = variableNames;
        this.setColorSame(defaultColor);
    }
    
    public BoxPlot (double[][] values, String[] variableNames) {
        valueList = new ArrayList();
        for (int i = 0; i < values.length; i++) {
            valueList.add(values[i]);
        }
        this.variableNames = variableNames;
        this.setColorSame(defaultColor);
    }
    
    public BoxPlot (int[][] values, String[] variableNames) {
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
    
    public BoxPlot (ArrayList<double[]> valueList, String[] variableNames) {
        this.valueList = valueList;
        this.variableNames = variableNames;
    }
    
   
    public void setColorRainbow () {
        this.colors = new String[variableNames.length];
        for (int i = 0; i < colors.length; i++) colors[i] = defaultColorArray[i];
    }
    
    public void setColorSame (String c) {
        this.colors = new String[variableNames.length];
        for (int i = 0; i < colors.length; i++) colors[i] = c;
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
    
    public void setXLabVertical () {
        this.ifXVertical = true;
    }
    
    @Override
    public void showGraph() {
        RCaller caller = new RCaller();
        caller.setRscriptExecutable(this.rPath);
        RCode rCode = new RCode();
        for (int i = 0; i < this.valueList.size(); i++) {
            rCode.addDoubleArray("x"+String.valueOf(i), this.valueList.get(i));
        }
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
        for (int i = 0; i < this.valueList.size(); i++) {
            rCode.addDoubleArray("x"+String.valueOf(i), this.valueList.get(i));
        }
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
        sb.append("boxplot(");
        for (int i = 0; i < this.valueList.size(); i++) {
            sb.append("x"+String.valueOf(i)).append(",");
        }
        if (this.ifXVertical) sb.append("las=2,");
        sb.append("col=c(");
        for (int i = 0; i < this.variableNames.length; i++) {
            sb.append("\"").append(colors[i]).append("\",");
        }
        sb.deleteCharAt(sb.length()-1);
        sb.append("),");
        sb.append("names=c(");
        for (int i = 0; i < this.variableNames.length; i++) {
            sb.append("\"").append(variableNames[i]).append("\",");
        }
        sb.deleteCharAt(sb.length()-1);
        sb.append("),");
        sb.append("main=\"").append(title).append("\",");
        sb.append("ylab=\"").append(yLab).append("\",");
        sb.append("cex.main=").append(titleSize).append(",");
        sb.append("cex.lab=").append(labSize).append(",");
        sb.append("cex.axis=").append(axisSize).append(",");
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
