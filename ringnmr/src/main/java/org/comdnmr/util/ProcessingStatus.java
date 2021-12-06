/*
 * CoMD/NMR Software : A Program for Analyzing NMR Dynamics Data
 * Copyright (C) 2018-2019 Bruce A Johnson
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
package org.comdnmr.util;

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
