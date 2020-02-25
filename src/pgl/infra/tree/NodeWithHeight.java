/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.infra.tree;

/**
 *
 * @author feilu
 */
public class NodeWithHeight implements Comparable<NodeWithHeight> {
    
    String name = null;
    
    double height = 0;
    
    public NodeWithHeight (String name, double height) {
        this.name = name;
        this.height = height;
    }
    
    public String getName () {
        return name;
    }
    
    public double getHeight () {
        return height;
    }
    
    public void setName (String name) {
        this.name = name;
    }
    
    public void setHeight (double height) {
        this.height = height;
    }
    
    @Override
    public int compareTo(NodeWithHeight o) {
        return name.compareTo(o.name);
//        if (this.height < o.height) return -1;
//        else if (this.height > o.height) return 1;
//        return 0;
    }
}
