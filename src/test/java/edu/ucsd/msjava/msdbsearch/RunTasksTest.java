package edu.ucsd.msjava.msdbsearch;

import org.junit.Test;

import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;

import static org.junit.Assert.*;

/**
 * ConcurrentMSGFPlus.runTasks has to report the failure that stopped a phase, not a failure of the clean-up after it.
 */
public class RunTasksTest {

    @Test
    public void runsEveryTask() throws Throwable {
        final AtomicInteger numTasksRun = new AtomicInteger();
        List<Runnable> tasks = new ArrayList<Runnable>();
        for (int i = 0; i < 20; i++) {
            tasks.add(new Runnable() {
                public void run() {
                    numTasksRun.incrementAndGet();
                }
            });
        }
        ConcurrentMSGFPlus.runTasks("Test", tasks, 4);
        assertEquals(20, numTasksRun.get());
    }

    @Test
    public void rethrowsTheExceptionOfATask() {
        final IllegalStateException failure = new IllegalStateException("task failure");
        List<Runnable> tasks = new ArrayList<Runnable>();
        tasks.add(new Runnable() {
            public void run() {
                throw failure;
            }
        });
        try {
            ConcurrentMSGFPlus.runTasks("Test", tasks, 1);
            fail("the task's exception was lost");
        } catch (Throwable e) {
            assertSame(failure, e);
        }
    }

    // The pool is shut down by the failure, so submitting the second task fails; the task's exception is the one to report
    @Test
    public void rethrowsTheExceptionOfATaskThatFailedWhileTasksWereSubmitted() {
        final IllegalStateException failure = new IllegalStateException("task failure");
        final CountDownLatch firstTaskStarted = new CountDownLatch(1);
        List<Runnable> tasks = new AbstractList<Runnable>() {
            public int size() {
                return 2;
            }

            public Runnable get(int index) {
                if (index == 0) {
                    return new Runnable() {
                        public void run() {
                            firstTaskStarted.countDown();
                            throw failure;
                        }
                    };
                }
                try {
                    firstTaskStarted.await();
                    Thread.sleep(200);
                } catch (InterruptedException e) {
                    throw new RuntimeException(e);
                }
                return new Runnable() {
                    public void run() {
                    }
                };
            }
        };
        try {
            ConcurrentMSGFPlus.runTasks("Test", tasks, 1);
            fail("the task's exception was lost");
        } catch (Throwable e) {
            assertSame(failure, e);
        }
    }

    // The interrupted driver stops its tasks, and a task that ends normally afterwards must not replace the interruption
    // with a failure to schedule its progress report
    @Test
    public void interruptingTheDriverStopsTheTasksAndIsReported() throws Exception {
        final CountDownLatch taskStarted = new CountDownLatch(1);
        final AtomicBoolean taskInterrupted = new AtomicBoolean(false);
        final AtomicReference<Throwable> thrown = new AtomicReference<Throwable>();
        final List<Runnable> tasks = Collections.<Runnable>singletonList(new Runnable() {
            public void run() {
                taskStarted.countDown();
                try {
                    new CountDownLatch(1).await();
                } catch (InterruptedException e) {
                    taskInterrupted.set(true);
                }
            }
        });
        Thread driver = new Thread(new Runnable() {
            public void run() {
                try {
                    ConcurrentMSGFPlus.runTasks("Test", tasks, 1);
                } catch (Throwable e) {
                    thrown.set(e);
                }
            }
        });
        driver.start();
        assertTrue(taskStarted.await(30, TimeUnit.SECONDS));
        driver.interrupt();
        driver.join(30000);
        assertFalse(driver.isAlive());
        assertTrue(taskInterrupted.get());
        assertTrue(String.valueOf(thrown.get()), thrown.get() instanceof InterruptedException);
    }

    // An interruption while the driver waits for its tasks to stop must not replace the failure that stopped them
    @Test
    public void interruptingTheCleanUpKeepsTheFirstFailure() throws Exception {
        final IllegalStateException failure = new IllegalStateException("task failure");
        final CountDownLatch stubbornTaskStarted = new CountDownLatch(1);
        // once by the pool when the other task fails, once by the driver's clean-up
        final CountDownLatch stubbornTaskInterruptedTwice = new CountDownLatch(2);
        final CountDownLatch releaseStubbornTask = new CountDownLatch(1);
        final AtomicReference<Throwable> thrown = new AtomicReference<Throwable>();
        final AtomicBoolean driverInterruptedAfterwards = new AtomicBoolean(false);
        final List<Runnable> tasks = new AbstractList<Runnable>() {
            public int size() {
                return 3;
            }

            public Runnable get(int index) {
                if (index == 0) {
                    return new Runnable() {
                        // ignores interruptions, so that the driver has to wait for it
                        public void run() {
                            stubbornTaskStarted.countDown();
                            while (true) {
                                try {
                                    releaseStubbornTask.await();
                                    return;
                                } catch (InterruptedException e) {
                                    stubbornTaskInterruptedTwice.countDown();
                                }
                            }
                        }
                    };
                }
                if (index == 1) {
                    return new Runnable() {
                        public void run() {
                            try {
                                stubbornTaskStarted.await();
                            } catch (InterruptedException e) {
                                throw new RuntimeException(e);
                            }
                            throw failure;
                        }
                    };
                }
                // submitted after the failure has shut the pool down, which sends the driver into its clean-up
                while (stubbornTaskInterruptedTwice.getCount() == 2)
                    Thread.yield();
                return new Runnable() {
                    public void run() {
                    }
                };
            }
        };
        Thread driver = new Thread(new Runnable() {
            public void run() {
                try {
                    ConcurrentMSGFPlus.runTasks("Test", tasks, 2);
                } catch (Throwable e) {
                    thrown.set(e);
                }
                driverInterruptedAfterwards.set(Thread.currentThread().isInterrupted());
            }
        });
        driver.start();
        assertTrue(stubbornTaskInterruptedTwice.await(30, TimeUnit.SECONDS));
        driver.interrupt();
        driver.join(30000);
        releaseStubbornTask.countDown();
        assertFalse(driver.isAlive());
        assertSame(failure, thrown.get());
        assertTrue(driverInterruptedAfterwards.get());
    }
}
