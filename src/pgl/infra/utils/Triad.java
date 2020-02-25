/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.infra.utils;

/**
 *
 * @author feilu
 * @param <T>
 * @param <S>
 * @param <U>
 */
public class Triad<T, S, U>  {
    T first = null;
    S second = null;
    U third = null;

    public Triad (T first, S second, U third) {
        this.first = first;
        this.second = second;
        this.third = third;
    }
    
    public T getFirstElement () {
        return first;
    }
    
    public S getSecondElement () {
        return second;
    }
    
    public U getThirdElement () {
        return third;
    }
}
