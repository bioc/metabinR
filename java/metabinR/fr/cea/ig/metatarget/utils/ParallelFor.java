/*
 *
 * MetaTarget ParallelFor
 *
 * Copyright (C) 2022 Anestis Gkanogiannis <anestis@gkanogiannis.com>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 *
 */
package fr.cea.ig.metatarget.utils;

import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;

public class ParallelFor {

	private static final int NUM_CORES = Runtime.getRuntime().availableProcessors();

	private static final ForkJoinPool fjPool = new ForkJoinPool(NUM_CORES, ForkJoinPool.defaultForkJoinWorkerThreadFactory, null, true);

	public static <T> void blockingFor(final Iterable<? extends T> elements, final Operation<T> operation) {
		blockingFor(2 * NUM_CORES, elements, operation);
	}

	public static <T> void blockingFor(int numThreads, final Iterable<? extends T> elements, final Operation<T> operation) {
		For(numThreads, new NamedThreadFactory("Parallel.For"), elements, operation, Integer.MAX_VALUE, TimeUnit.DAYS);
	}

	public static <T> void For(final Iterable<? extends T> elements, final Operation<T> operation) {
		For(2 * NUM_CORES, elements, operation);
	}

	public static <T> void For(int numThreads, final Iterable<? extends T> elements, final Operation<T> operation) {
		For(numThreads, new NamedThreadFactory("Parallel.For"), elements, operation, null, null);
	}

	public static <S extends T, T> void For(int numThreads, NamedThreadFactory threadFactory, final Iterable<S> elements, final Operation<T> operation, Integer wait, TimeUnit waitUnit) {
		ThreadPoolExecutor threadPoolExecutor = new ThreadPoolExecutor(numThreads, numThreads, 0L, TimeUnit.MILLISECONDS,new LinkedBlockingQueue<Runnable>());
		final ThreadSafeIterator<S> itr = new ThreadSafeIterator<S>(elements.iterator());

		for (int i = 0; i < threadPoolExecutor.getMaximumPoolSize(); i++) {
			threadPoolExecutor.submit(new Callable<Void>() {
				@Override
				public Void call() {
					T element;
					while ((element = itr.next()) != null) {
						try {
							operation.perform(element);
						} catch (Exception e) {
							Logger.getLogger(ParallelFor.class.getName()).log(Level.SEVERE, "Exception during execution of parallel task", e);
						}
					}
					return null;
				}
			});
		}

		threadPoolExecutor.shutdown();

		if (wait != null) {
			try {
				threadPoolExecutor.awaitTermination(wait, waitUnit);
			} catch (InterruptedException e) {
				throw new IllegalStateException(e);
			}
		}
	}

	private static class ThreadSafeIterator<T> {

		private final Iterator<T> itr;

		public ThreadSafeIterator(Iterator<T> itr) {
			this.itr = itr;
		}

		public synchronized T next() {
			return itr.hasNext() ? itr.next() : null;
		}
	}

	public static <T> void ForFJ(final Iterable<T> elements, final Operation<T> operation) {
		fjPool.invokeAll(createCallables(elements, operation));
	}

	public static <T> Collection<Callable<Void>> createCallables(final Iterable<T> elements, final Operation<T> operation) {
		List<Callable<Void>> callables = new LinkedList<Callable<Void>>();
		for (final T elem : elements) {
			callables.add(new Callable<Void>() {

				@Override
				public Void call() {
					operation.perform(elem);
					return null;
				}
			});
		}

		return callables;
	}

	public static interface Operation<T> {
		public void perform(T pParameter);
	}
	
}

/*
 	public class Parallel {
    private static final int NUM_CORES = Runtime.getRuntime().availableProcessors();

    private static final ExecutorService forPool = Executors.newFixedThreadPool(NUM_CORES * 2, new NamedThreadFactory("Parallel.For"));

    public static <T> void For(final Iterable<T> elements, final Operation<T> operation) {
        try {
            // invokeAll blocks for us until all submitted tasks in the call complete
            forPool.invokeAll(createCallables(elements, operation));
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }

    public static <T> Collection<Callable<Void>> createCallables(final Iterable<T> elements, final Operation<T> operation) {
        List<Callable<Void>> callables = new LinkedList<Callable<Void>>();
        for (final T elem : elements) {
            callables.add(new Callable<Void>() {
                @Override
                public Void call() {
                    operation.perform(elem);
                    return null;
                }
            });
        }

        return callables;
    }

    public static interface Operation<T> {
        public void perform(T pParameter);
    }
}
 */

