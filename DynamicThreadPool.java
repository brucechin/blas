import com.sun.management.OperatingSystemMXBean;
import java.lang.management.ManagementFactory;
import java.util.Random;
import java.util.Timer;
import java.util.TimerTask;
import java.util.concurrent.*;
import org.apache.log4j.Logger;
import org.apache.log4j.LogManager;

/**
 * @CLassName DynamicThreadPool
 * @Description modify number of threads in the thread pool according to the utilization of host machine CPU
 * @Author lianke.qin@gmail.com
 * @Date 2019/7/25 2:44 PM
 * @Version 1.0
 **/
public class DynamicThreadPool<Job extends Runnable>{
    private int numProcessors;
    private double procUtil;//utilization percentage of all processors
    private int minThreads;
    private int maxThreads;
    private int curThreads;//current number of threads
    private int updatePeriod;//run updateWorkers() once per period
    public ThreadPoolExecutor threadPool;
    OperatingSystemMXBean bean;
    Timer timer;
    private static Logger logger = Logger.getLogger(DynamicThreadPool.class);
    /*
    *
    * @Param minthread : min number of threads in this DynamicThreadPool
    * @Param maxthread : max number of threads in this DynamicThreadPool
    * @Param period : run updateWorkers() every period milliseconds
    * @Param queue : used to init the waiting queue in ThreadPool
    *
    * */
    public DynamicThreadPool(int minthread, int maxthread, int period, BlockingQueue queue){
        bean = (OperatingSystemMXBean) ManagementFactory.getOperatingSystemMXBean();
        numProcessors = bean.getAvailableProcessors();
        procUtil = bean.getSystemCpuLoad();
        minThreads = minthread;
        maxThreads = maxthread;
        curThreads = minThreads;
        Logger logger = Logger.getLogger(DynamicThreadPool.class);
        updatePeriod = period;
        timer = new Timer();
        timer.schedule(new ThreadTask(), 100, updatePeriod);
        threadPool = new ThreadPoolExecutor(minThreads, maxThreads, 10000, TimeUnit.MILLISECONDS, queue);
    }

    public void setUpdatePeriod(int period){
        updatePeriod = period;
    }

    public void printInfo(){
        logger.info("CPU Utilization : " + procUtil);
        //System.out.println("CPU Utilization : " + procUtil);
        System.out.println("Current size of thread pool : " + curThreads);
    }

    public void updateWorkers(){

        procUtil = bean.getSystemCpuLoad();
        if(curThreads != minThreads + (int)((maxThreads - minThreads) * (1 - procUtil))) {
            //reduce number of setCorePoolSize calls
            curThreads = minThreads + (int) ((maxThreads - minThreads) * (1 - procUtil));
            threadPool.setCorePoolSize(curThreads);
        }
    }

    public void execute(Job j){
        threadPool.execute(j);
    }

    public long getCompletedTaskCount(){
        return threadPool.getCompletedTaskCount();
    }

    public long getTaskCount(){
        return threadPool.getTaskCount();
    }

    public int getPoolSize(){
        return threadPool.getPoolSize();
    }

    public void shutdown(){
        //stop receiving new tasks, but execute the tasks already in queue
        timer.cancel();
        threadPool.shutdown();
        System.out.println("Dynamic thread pool is shutting down, no more new tasks accepted.");
        try {//wait for tasks in queue to complete
            threadPool.awaitTermination(Long.MAX_VALUE, TimeUnit.MINUTES);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        System.out.println("Dynamic thread pool terminated, all tasks are executed");

    }

    public boolean isShutdown(){
        return threadPool.isShutdown();
    }

    public void shutdownNow(){
        //reject new incoming tasks, also stop executing tasks in the queue
        timer.cancel();
        threadPool.shutdownNow();
        System.out.println("Dynamic thread pool shutdowns forcefully ");
    }

    public int getMaxThreads(){
        return maxThreads;
    }

    public int getMinThreads(){
        return minThreads;
    }

    public void setMinThreads(int minthreads){
        minThreads = minthreads;
    }

    public void setMaxThreads(int maxthreads){
        maxThreads = maxthreads;
    }


    public class ThreadTask extends TimerTask{
        public ThreadTask(){
        }

        @Override
        public void run(){
            updateWorkers();
            printInfo();
        }
    }

    public static void main(String[] args){
        BlockingQueue queue = new ArrayBlockingQueue<>(500);
        DynamicThreadPool test = new DynamicThreadPool(2,5, 2000, queue);

        for(int i = 0; i < 1000; i++){
            if(i == 500) test.shutdownNow();
            logger.info(" Task " + test.threadPool.getTaskCount() + " inserted to the queue");
            try{
                Thread.sleep((int)(Math.random() * 20) + 20);
                test.execute(new RealTask(i));
            }catch (Exception e){
                e.printStackTrace();
            }

        }
        try{
            Thread.sleep(10000);
        }catch (Exception e){
            e.printStackTrace();
        }finally {
            test.shutdown();
        }

    }

}
