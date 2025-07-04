import time
import pyrbd_plusplus as pyrbd

def benchmark_optimization():
    # 准备大规模测试数据
    node_pairs = [(i, i+1) for i in range(100)]
    probabilities = {i: 0.95 for i in range(200)}
    pathsets_list = [[[i, i+1, i+2] for i in range(j, j+50)] for j in range(100)]
    
    # 测试串行版本
    start = time.time()
    result1 = pyrbd.sdp.eval_avail_topo(node_pairs, probabilities, pathsets_list)
    serial_time = time.time() - start
    
    # 测试并行版本
    start = time.time()
    result2 = pyrbd.sdp.eval_avail_topo_parallel(node_pairs, probabilities, pathsets_list)
    parallel_time = time.time() - start
    
    speedup = serial_time / parallel_time
    print(f"Serial time: {serial_time:.3f}s")
    print(f"Parallel time: {parallel_time:.3f}s")
    print(f"Speedup: {speedup:.2f}x")

if __name__ == "__main__":
    benchmark_optimization()