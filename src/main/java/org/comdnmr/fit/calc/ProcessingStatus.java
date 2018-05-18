package org.comdnmr.fit.calc;

/**
 *
 * @author Bruce Johnson
 */
public class ProcessingStatus {

    /**
     * @return the status
     */
    public String getStatus() {
        return status;
    }

    /**
     * @return the ok
     */
    public boolean isOk() {
        return ok;
    }

    /**
     * @return the throwable
     */
    public Throwable getThrowable() {
        return throwable;
    }

    private final String status;
    private final boolean ok;
    private final Throwable throwable;

    public ProcessingStatus(String status) {
        this(status, true, null);
    }

    public ProcessingStatus(String status, boolean ok) {
        this(status, ok, null);
    }

    public ProcessingStatus(String status, boolean ok, Throwable throwable) {
        this.status = status;
        this.ok = ok;
        this.throwable = throwable;
    }

}
