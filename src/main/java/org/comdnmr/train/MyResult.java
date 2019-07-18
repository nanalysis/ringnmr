package org.comdnmr.fit.train;

/**
 * Purpose of this class is to allow users to retrieve an object with references
 * to at most 3 other objects of different types. 
 *
 * @author Teddy Colon
 *
 * @param <F> First object to return of type any type
 * @param <S> Second object to return of type any type
 * @param <T> Third object to return of type any type
 */
public final class MyResult<F, S, T> {

    private final F first;
    private final S second;
    private final T third;
    
    
    public MyResult(F first, S second, T third) {
        this.first = first;
        this.second = second;
        this.third = third;
    }
    
    public MyResult cloning() {
        return new MyResult(first,second,third);
    }

    public F getFirst() {
        return first;
    }

    public S getSecond() {
        return second;
    }

    public T getThird() {
        return third;
    }

}
