/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.infra.tree;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.List;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.TreeNode;

/**
 *
 * @author feilu
 */
public class Newick {
    
    DefaultMutableTreeNode root = new DefaultMutableTreeNode(new NodeWithHeight("root", 0));
    List<NodeWithHeight> leafList = new ArrayList<>();
    List<String> taxaList = new ArrayList();
    HashMap<String, DefaultMutableTreeNode> taxaNodeMap = new HashMap<>();
    
    public Newick (String nwkS) {
        this.readNwk(nwkS, root);
        this.buildTaxaListAndMap();
    }
    
    /**
     * 
     * @param n
     * @return 
     */
    public List<String> selectTaxaWithMaxDiversity (int n) {
        while (root.getLeafCount() > n) {
            DefaultMutableTreeNode minDmt = this.getLeafNodeWithMinHeight();
            DefaultMutableTreeNode parentDmt = (DefaultMutableTreeNode)minDmt.getParent();
            minDmt.removeFromParent();
            if (parentDmt.getChildCount() == 1) {
                NodeWithHeight parentNwh = (NodeWithHeight)parentDmt.getUserObject();
                DefaultMutableTreeNode childDmt = (DefaultMutableTreeNode)parentDmt.getFirstChild();
                NodeWithHeight childNwh = (NodeWithHeight)childDmt.getUserObject();                
                double newParentHeight = parentNwh.height+childNwh.height;
                DefaultMutableTreeNode newParent = new DefaultMutableTreeNode(new NodeWithHeight(childNwh.name, newParentHeight));
                if (!childDmt.isLeaf()) {
                    newParent.add(childDmt);
                }
                DefaultMutableTreeNode grandDmt = (DefaultMutableTreeNode)parentDmt.getParent();
                parentDmt.removeFromParent();
                if (grandDmt != null) grandDmt.add(newParent); 
                else root = newParent;
            }
        }
        leafList.clear();
        taxaList.clear();
        taxaNodeMap.clear();
        this.buildTaxaListAndMap();
        return this.taxaList;
    }
    
    public DefaultMutableTreeNode getLeafNodeWithMinHeight () {
        DefaultMutableTreeNode currentDmt = null;
        double currentV = Double.MAX_VALUE;
        Enumeration e = root.breadthFirstEnumeration();
        while (e.hasMoreElements()) {
            DefaultMutableTreeNode d = (DefaultMutableTreeNode)e.nextElement();
            if (!d.isLeaf()) continue;
            NodeWithHeight ee = (NodeWithHeight)d.getUserObject();
            if (ee.getHeight() < currentV) {
                currentV = ee.getHeight();
                currentDmt = d;
            }
        }
        return currentDmt;
    }
    
    public String getTaxonWithMinHeight () {
        String currentT = null;
        double currentV = Double.MAX_VALUE;
        for (int i = 0; i < taxaList.size(); i++) {
            DefaultMutableTreeNode dmt = this.getTaxonNode(taxaList.get(i));
            NodeWithHeight nwh = (NodeWithHeight)dmt.getUserObject();
            if (nwh.getHeight() < currentV) {
                currentT = nwh.getName();
                currentV = nwh.getHeight();
            }
        }
        return currentT;
    }
    
    public double getHeight (String taxon) {
        DefaultMutableTreeNode dmt = this.getTaxonNode(taxon);
        NodeWithHeight nwh = (NodeWithHeight)dmt.getUserObject();
        return nwh.height;
    }
    
    public double getMaxHeightToRootAcrossTaxa () {
        double max = -1;
        for (int i = 0; i < taxaList.size(); i++) {
            double cHeight = this.getHeightToRoot(taxaList.get(i));
            if (cHeight > max) max = cHeight;
        }
        return max;
    }
    
    public double getHeightToRoot (String taxon) {
        double[] height = this.getPathToRootOfHeight(taxon);
        double sum = 0;
        for (int i = 0; i < height.length; i++) {
            sum+=height[i];
        }
        return sum;
    }
    
    public double getHeightToRoot (DefaultMutableTreeNode node) {
        double[] height = this.getPathToRootOfHeight(node);
        double sum = 0;
        for (int i = 0; i < height.length; i++) {
            sum+=height[i];
        }
        return sum;
    }
    
    public double[] getPathToRootOfHeight (DefaultMutableTreeNode node) {
        NodeWithHeight[] nhPath = this.getPathToRoot(node);
        double[] heightPath = new double[nhPath.length];
        for (int i = 0; i < heightPath.length; i++) {
            heightPath[i] = nhPath[i].getHeight();
        }
        return heightPath;
    }
    
    public String[] getPathToRootOfTaxa (DefaultMutableTreeNode node) {
        NodeWithHeight[] nhPath = this.getPathToRoot(node);
        String[] taxaPath = new String[nhPath.length];
        for (int i = 0; i < taxaPath.length; i++) {
            taxaPath[i] = nhPath[i].getName();
            System.out.println(taxaPath[i]);
        }
        return taxaPath;
    }
    
    public double[] getPathToRootOfHeight (String taxon) {
        NodeWithHeight[] nhPath = this.getPathToRoot(taxon);
        double[] heightPath = new double[nhPath.length];
        for (int i = 0; i < heightPath.length; i++) {
            heightPath[i] = nhPath[i].getHeight();
        }
        return heightPath;
    }
    
    public String[] getPathToRootOfTaxa (String taxon) {
        NodeWithHeight[] nhPath = this.getPathToRoot(taxon);
        String[] taxaPath = new String[nhPath.length];
        for (int i = 0; i < taxaPath.length; i++) {
            taxaPath[i] = nhPath[i].getName();
            System.out.println(taxaPath[i]);
        }
        return taxaPath;
    }
    
    private NodeWithHeight[] getPathToRoot (String taxon) {
        DefaultMutableTreeNode tnode = this.getTaxonNode(taxon);
        TreeNode[] path = tnode.getPath();
        NodeWithHeight[] nhPath = new NodeWithHeight[path.length];
        for (int i = 0; i < path.length; i++) {
            DefaultMutableTreeNode no = (DefaultMutableTreeNode)path[i];
            NodeWithHeight nh = (NodeWithHeight)no.getUserObject();
            nhPath[i] = nh;
        }
        return nhPath;
    }
    
    private NodeWithHeight[] getPathToRoot (DefaultMutableTreeNode node) {
        TreeNode[] path = node.getPath();
        NodeWithHeight[] nhPath = new NodeWithHeight[path.length];
        for (int i = 0; i < path.length; i++) {
            DefaultMutableTreeNode no = (DefaultMutableTreeNode)path[i];
            NodeWithHeight nh = (NodeWithHeight)no.getUserObject();
            nhPath[i] = nh;
        }
        return nhPath;
    }
    
    //need to adjust corresponding fields
    public boolean deleteTaxonNode (String taxon) {
        DefaultMutableTreeNode  node = this.getTaxonNode(taxon);
        if (node == null) return false;
        else {
            node.removeFromParent();
        }
        return true;
    }
    
    public DefaultMutableTreeNode getTaxonNode (String taxon) {
        return taxaNodeMap.get(taxon);
    }
    
    public double getDistanceBetweenTaxa (String taxon1, String taxon2) {
        DefaultMutableTreeNode sNode = this.getSharedAncesterNode(taxon1, taxon2);
        double dis = this.getHeightToRoot(taxon1) + this.getHeightToRoot(taxon2)-2*this.getHeightToRoot(sNode);
        return dis;
    }
    
    public String getSharedAncesterTaxon (String taxon1, String taxon2) {
        DefaultMutableTreeNode node = this.getSharedAncesterNode(taxon1, taxon2);
        return Newick.getTaxonName(node);
    }
    
    public DefaultMutableTreeNode getSharedAncesterNode (String taxon1, String taxon2) {
        DefaultMutableTreeNode tnode = this.getTaxonNode(taxon1);
        TreeNode t = tnode.getSharedAncestor(this.getTaxonNode(taxon2));
        DefaultMutableTreeNode sNode= (DefaultMutableTreeNode)t;
        return sNode;
    }
    
    public static String getTaxonName (DefaultMutableTreeNode node) {
        return ((NodeWithHeight)node.getUserObject()).getName();
    }
    
    public static double getHeight (DefaultMutableTreeNode node) {
        return ((NodeWithHeight)node.getUserObject()).getHeight();
    }
    
    public int getTaxaCount () {
        return taxaList.size();
    }
    
    private void buildTaxaListAndMap () {
        Enumeration e = root.breadthFirstEnumeration();
        while (e.hasMoreElements()) {
            DefaultMutableTreeNode d = (DefaultMutableTreeNode)e.nextElement();
            if (!d.isLeaf()) continue;
            NodeWithHeight ee = (NodeWithHeight)d.getUserObject();
            leafList.add(ee);
            taxaNodeMap.put(ee.getName(), d);
        }
        Collections.sort(leafList);
        for (int i = 0; i < leafList.size(); i++) {
            taxaList.add(leafList.get(i).getName());
        }
    }
    
    private void readNwk (String nwkS, DefaultMutableTreeNode parent) {
        if (nwkS.endsWith(";")) nwkS = nwkS.replaceFirst(";", "");
        if (nwkS.startsWith("(") && nwkS.endsWith(")")) {
            nwkS = nwkS.substring(1, nwkS.length()-1);
        }
        int cnt = 0;
        int currentIndex = 0;
        for (int i = 0; i < nwkS.length(); i++) {
            if (nwkS.charAt(i) == '(') cnt++;
            else if (nwkS.charAt(i) == ')') cnt--;
            else if (nwkS.charAt(i) == ',') {
                if (cnt == 0) {
                    String currentNwkS = nwkS.substring(currentIndex, i);
                    currentIndex = i+1;
                    //System.out.println(currentNwkS);                    
                    int cIndex = currentNwkS.length()-1;
                    for (int j = 0; j < currentNwkS.length(); j++) {
                        if (currentNwkS.charAt(cIndex) == ':') break;
                        cIndex--;
                    }
                    String name = currentNwkS.substring(0, cIndex);
                    double height = Double.parseDouble(currentNwkS.substring(cIndex+1, currentNwkS.length()));
                    DefaultMutableTreeNode child = new DefaultMutableTreeNode(new NodeWithHeight(name, height));
                    parent.add(child);
                    if (name.contains(")")){
                        this.readNwk(name, child);
                    }                                                
                }
            }
            if (i == nwkS.length()-1) {
                String currentNwkS = nwkS.substring(currentIndex, i+1);
                //System.out.println(currentNwkS);
                int cIndex = currentNwkS.length()-1;
                for (int j = 0; j < currentNwkS.length(); j++) {
                    if (currentNwkS.charAt(cIndex) == ':') break;
                    cIndex--;
                }
                String name = currentNwkS.substring(0, cIndex);
                double height = Double.parseDouble(currentNwkS.substring(cIndex+1, currentNwkS.length()));
                DefaultMutableTreeNode child = new DefaultMutableTreeNode(new NodeWithHeight(name, height));
                parent.add(child);
                if (name.contains(")")){
                    this.readNwk(name, child);
                }                     
            }
        }
    }
}
