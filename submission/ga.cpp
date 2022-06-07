#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cerrno>
#include <cstring>

#include <vector>
#include <queue>
#include <tuple>
#include <algorithm>

#define MAX_V 5000
#define MAX_E 40000
#define MOD 1000000007
#define SPARE_TIME 1

#define MAX_POPULATION 32768
#define NUM_CHILDREN 1024
#define CUTTING_POINT 2

double starts_at;

int V, E;
int edges[MAX_E][3];
std::vector<int> vertices[MAX_V];
std::vector<std::pair<int, int>> infos[MAX_V];
std::priority_queue<std::pair<int, int>> Q;

int hash_const;  // (2^V - 1) % MOD
int renumber_cnt;
uint8_t visits[MAX_V];
int degrees[MAX_V], renumbers[MAX_V], real_numbers[MAX_V];

int get_complement(int hash) {
	int complement = hash_const - hash;
	if (complement < 0) complement += MOD;
	return complement;
}

double get_time() {
	struct timespec ts;
	// get wall-clock time
	if (clock_gettime(CLOCK_REALTIME, &ts)) exit(errno);
	return ts.tv_sec + ts.tv_nsec * 1e-9;
}

class chromosome {
public:
	uint8_t genes[MAX_V];
	int score = INT32_MAX;

	chromosome(bool initialize = false) {
		if (initialize)
			for (int i = 0; i < V; i += 31) {
				int rand_num = rand();  // 31-bits
				int end = std::max(i + 31, V);
				for (int j = i; j < end; j++) {
					genes[j] = rand_num % 2;
					rand_num >>= 1;
				}
			}
	}

	chromosome(chromosome *other) {
		std::memcpy(genes, other->genes, V * sizeof(uint8_t));
	}

	void get_interval(chromosome *other, int left, int right, bool flip) {
		std::memcpy(genes + left, other->genes + left, (right - left) * sizeof(uint8_t));
		if (flip)
			for (int i = left; i < right; i++)
				genes[i] = 1 - genes[i];
	}

	chromosome *crossover(chromosome *other) {
		int cp[CUTTING_POINT + 2];
		cp[0] = 0;
		cp[CUTTING_POINT + 1] = V;
		for (int i = 1; i <= CUTTING_POINT; i++)
			cp[i] = rand() % V;
		// std::sort(cp + 1, cp + CUTTING_POINT + 1);
		if (cp[1] > cp[2]) {
			int temp = cp[1];
			cp[1] = cp[2];
			cp[2] = temp;
		}
		// create empty chromosome and copy intervals from this and others
		chromosome *child = new chromosome();
		bool flip = rand() % 2;
		for (int i = 0; i <= CUTTING_POINT; i++)
			child->get_interval((i % 2) ? other : this, cp[i], cp[i + 1], (i % 2) && flip);
		return child;
	}

	chromosome *mutation(bool create) {
		int idx = rand() % V;
		if (create) {
			chromosome *child = new chromosome(this);
			child->genes[idx] = 1 - child->genes[idx];
			return child;
		} else {
			genes[idx] = 1 - genes[idx];
			return this;
		}
	}

	chromosome *local_opt() {
		std::memset(degrees, 0, V * sizeof(int));
		score = 0;
		for (int i = 0; i < E; i++) {
			if (genes[edges[i][0]] != genes[edges[i][1]]) {
				score += edges[i][2];
				degrees[edges[i][0]] -= edges[i][2];
				degrees[edges[i][1]] -= edges[i][2];
			} else {
				degrees[edges[i][0]] += edges[i][2];
				degrees[edges[i][1]] += edges[i][2];
			}
		}
		for (int i = 0; i < V; i++)
			if (degrees[i] > 0)
				Q.emplace(degrees[i], i);
		while (!Q.empty()) {
			auto [diff, u] = Q.top();
			Q.pop();
			if (diff != degrees[u])
				continue;
			score += diff;
			for (auto [v, w] : infos[u]) {
				if (genes[u] != genes[v]) {
					degrees[u] += 2 * w;
					degrees[v] += 2 * w;
				} else {
					degrees[u] -= 2 * w;
					degrees[v] -= 2 * w;
				}
				if (degrees[v] > 0)
					Q.emplace(degrees[v], v);
			}
			genes[u] = 1 - genes[u];
			if (degrees[u] > 0)
				Q.emplace(degrees[u], u);
		}
		return this;
	}

	int hash() const {
		int result = 0;
		for (int i = 0; i < V; i++)
			result = ((result << 1) | genes[i]) % MOD;
		return std::min(result, get_complement(result));
	}

	int evaluate() {
		if (score == INT32_MAX) {
			score = 0;
			for (int i = 0; i < E; i++)
				if (genes[edges[i][0]] != genes[edges[i][1]])
					score += edges[i][2];
		}
		return score;
	}
};

class evaluation {
public:
	int score;
	int hash;
	chromosome *chr;

	evaluation() = default;
	evaluation(chromosome *chr_) : chr(chr_) {
		score = chr->evaluate();
		hash = chr->hash();
	}

	bool operator<(const evaluation &other) const {
		return score > other.score ||
		       (score == other.score && hash < other.hash);
	}
};

class population {
public:
	int num_chrs;
	chromosome *chrs[MAX_POPULATION], *children[NUM_CHILDREN];
	evaluation evals[MAX_POPULATION + NUM_CHILDREN], temp[MAX_POPULATION + NUM_CHILDREN];

	population() {
		for (int i = 0; i < MAX_POPULATION; i++) {
			chrs[i] = new chromosome(true);
			evals[i] = evaluation(chrs[i]);
		}
		std::sort(evals, evals + MAX_POPULATION);
		num_chrs = 1;
		for (int i = 1; i < MAX_POPULATION; i++) {
			if (evals[i].score == evals[i - 1].score &&
			    evals[i].hash == evals[i - 1].hash)
				delete evals[i].chr;
			else
				evals[num_chrs++] = evals[i];
		}
		for (int i = 0; i < num_chrs; i++)
			chrs[i] = evals[i].chr;
	}

	void replace() {
		for (int i = 0; i < NUM_CHILDREN; i++)
			evals[num_chrs + i] = evaluation(children[i]);
		int total = num_chrs + NUM_CHILDREN;
		std::sort(evals + num_chrs, evals + total);
		int p = 0, q = num_chrs, r = 0;
		while (p < num_chrs || q < total) {
			if (p == num_chrs)
				temp[r] = evals[q++];
			else if (q == total)
				temp[r] = evals[p++];
			else {
				if (evals[p] < evals[q])
					temp[r] = evals[p++];
				else
					temp[r] = evals[q++];
			}
			r++;
		}
		evals[0] = temp[0];
		chrs[0] = evals[0].chr;
		num_chrs = 1;
		for (int i = 1; i < total; i++) {
			if (temp[i].score == temp[i - 1].score &&
			    temp[i].hash == temp[i - 1].hash)
				delete temp[i].chr;
			else {
				if (num_chrs < MAX_POPULATION) {
					evals[num_chrs] = temp[i];
					chrs[num_chrs] = temp[i].chr;
					num_chrs++;
				} else {
					delete temp[i].chr;
				}
			}
		}
	}
} *group;

void get_input() {
	// get input
	int u, v, w;
	if (scanf("%d %d", &V, &E) != 2) exit(errno);
	for (int i = 0; i < E; i++) {
		if (scanf("%d %d %d", &u, &v, &w) != 3) exit(errno);
		// change from 1-base to 0-base
		u--;
		v--;
		edges[i][0] = u;
		edges[i][1] = v;
		edges[i][2] = w;
		vertices[u].push_back(v);
		vertices[v].push_back(u);
	}

	// set hash constant (2^V - 1)
	hash_const = 1;
	for (int i = 0; i < V; i++)
		hash_const = (hash_const << 1) % MOD;
	hash_const = (hash_const + MOD - 1) % MOD;
}

void dfs(int u) {
	if (visits[u]) return;
	visits[u] = 1;
	renumbers[u] = renumber_cnt;
	real_numbers[renumber_cnt] = u;
	renumber_cnt++;
	for (int v : vertices[u])
		dfs(v);
}

void renumber() {
	Q.emplace(0, rand() % V);
	while (!Q.empty()) {
		int u = Q.top().second;
		Q.pop();
		if (visits[u]) continue;
		visits[u] = 1;
		renumbers[u] = renumber_cnt++;
		for (int v : vertices[u]) {
			if (visits[v]) continue;
			degrees[v]++;
			Q.emplace(degrees[v], v);
		}
	}
	for (int i = 0; i < V; i++)
		std::sort(vertices[i].begin(), vertices[i].end(),
			[&](int x, int y) {
				return renumbers[x] < renumbers[y];
			}
		);
	// concerns: some vertices cannot be visited
	// but seeming there is no such vertex, all seem to be connected
	std::memset(visits, 0, sizeof(visits));
	renumber_cnt = 0;
	for (int i = 0; i < V; i++)
		dfs(i);
	for (int i = 0; i < E; i++) {
		edges[i][0] = renumbers[edges[i][0]];
		edges[i][1] = renumbers[edges[i][1]];
		infos[edges[i][0]].emplace_back(edges[i][1], edges[i][2]);
		infos[edges[i][1]].emplace_back(edges[i][0], edges[i][2]);
	}
}

void try_GA() {
	group = new population();
	do {
		for (int i = 0; i < NUM_CHILDREN; i++) {
			int x = rand() % group->num_chrs;
			int y = rand() % group->num_chrs;
			group->children[i] = group->chrs[x]->crossover(group->chrs[y])  // always create new chromosome
											   ->mutation(false)  // create new chromosome when flag is given
											   ->local_opt();  // do not create new chromosome
		}
		group->replace();
	} while (get_time() - starts_at + SPARE_TIME < V / 6.0);
}

void print_output() {
	chromosome *best = group->evals[0].chr;

	for (int i = 0; i < V; i++)
		visits[real_numbers[i]] = best->genes[i];
	for (int i = 0; i < V; i++)
		if (visits[i])
			printf("%d ", i + 1);
	printf("\n");
}

int main() {
	// get start time
	starts_at = get_time();

	// srand, rand is fast, we do not need true-randomness
	srand(time(NULL));

	get_input();
	renumber();
	try_GA();
	print_output();
}
