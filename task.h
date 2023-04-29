
int task_create(void*(*task)(void *),int i, int period, int drel, int prio);

int get_task_index(void* arg);

void set_activation(int i);

int deadline_miss(int i);

void wait_for_activation(int i);

int task_delete(int i);

