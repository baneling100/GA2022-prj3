#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cerrno>

#include <vector>
#include <queue>
#include <algorithm>

#define MAX_V 5000
#define MAX_L 625  // L == (V + 7) / 8
#define MAX_E 40000
#define MOD 1000000007
#define SPARE_TIME 1

#define MAX_POPULATION 10240
#define NUM_CHILDREN 256
#define CUTTING_POINT 2

typedef std::pair<int, int> ipair;

double starts_at;

int V, L, E;
int edges[MAX_E][3];
std::vector<int> vertices[MAX_V];

int hash_const;  // (2^V - 1) % MOD

int renumber_cnt;
bool visits[MAX_V], answers[MAX_V];
int prenumbers[MAX_V], degrees[MAX_V], renumbers[MAX_V], real_numbers[MAX_V];
int renumbered_edges[MAX_E][3];

int get_complement(int hash) {
	int complement = hash_const - hash;
	if (complement < 0) complement += MOD;
	return complement;
}

class chromosome {
public:
	unsigned char genes[MAX_L];

	chromosome(bool is_random = false) {
		if (is_random) {
			for (int i = 0; i < L; i++)
				genes[i] = rand() % 256;
			// erase tail, if not, hash value can be wrong
			genes[L - 1] &= (1 << ((V - 1) % 8 + 1)) - 1;
		}
	}

	chromosome(chromosome *other) {
		for (int i = 0; i < L; i++)
			genes[i] = other->genes[i];
	}

	// copy interval [left, (right + 7) / 8 * 8) from other to this
	void get_interval(chromosome *other, int left, int right) {
		int pos1 = left / 8, pos2 = left % 8;
		genes[pos1] = (genes[pos1] & ((1 << pos2) - 1))
			    | (other->genes[pos1] & (255 << pos2));
		for (int i = pos1 + 1; i <= (right - 1) / 8; i++)
			genes[i] = other->genes[i];
	}

	chromosome *crossover(chromosome *other) {
		int cp[CUTTING_POINT + 2];
		cp[0] = 0;
		cp[CUTTING_POINT + 1] = V;
		for (int i = 1; i <= CUTTING_POINT; i++)
			cp[i] = rand() % V;
		std::sort(cp + 1, cp + CUTTING_POINT + 1);
		// create empty chromosome and copy intervals from this and others
		chromosome *child = new chromosome();
		for (int i = 0; i <= CUTTING_POINT; i++)
			child->get_interval((i % 2) ? other : this, cp[i], cp[i + 1]);
		return child;
	}

	chromosome *mutation() {
		chromosome *child = new chromosome(this);
		int idx = rand() % V;
		child->genes[idx / 8] ^= 1 << (idx % 8);
		return child;
	}

	int hash() const {
		long long result = 0LL;
		for (int i = 0; i < L; i++)
			result = (256LL * result + genes[i]) % MOD;
		return std::min((int)result, get_complement((int)result));
	}

	int evaluate() const {
		int score = 0;
		for (int i = 0; i < E; i++) {
			int u = renumbered_edges[i][0];
			int v = renumbered_edges[i][1];
			if ((genes[u / 8] & (1 << (u % 8))) !=
			    (genes[v / 8] & (1 << (v % 8))))
			    score += renumbered_edges[i][2];
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

	bool operator<(const evaluation &other) {
		return score > other.score ||
		       (score == other.score && hash < other.hash);
	}
};

class population {
public:
	int num_chrs;
	chromosome *chrs[MAX_POPULATION], *children[NUM_CHILDREN];
	evaluation evals[MAX_POPULATION + NUM_CHILDREN];

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
		std::sort(evals, evals + total);
		num_chrs = 1;
		for (int i = 1; i < total; i++) {
			if (evals[i].score == evals[i - 1].score &&
			    evals[i].hash == evals[i - 1].hash)
				delete evals[i].chr;
			else
				evals[num_chrs++] = evals[i];
		}
		total = num_chrs;
		num_chrs = std::min(num_chrs, MAX_POPULATION);
		for (int i = 0; i < num_chrs; i++)
			chrs[i] = evals[i].chr;
		for (int i = num_chrs; i < total; i++)
			delete evals[i].chr;
	}
} group;

double get_time() {
	struct timespec ts;
	// get wall-clock time
	if (clock_gettime(CLOCK_REALTIME, &ts)) exit(errno);
	return ts.tv_sec + ts.tv_nsec * 1e-9;
}

void get_input() {
	// get input
	int u, v, w;
	if (scanf("%d %d", &V, &E) != 2) exit(errno);
	L = (V + 7) / 8;
	for (int i = 0; i < E; i++) {
		if (scanf("%d %d %d", &u, &v, &w) != 3) exit(errno);
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
	for (int i = 0; i < V; i++) {
		hash_const *= 2;
		hash_const %= MOD;
	}
	hash_const += MOD - 1;
	hash_const %= MOD;
}

void dfs(int u) {
	if (visits[u]) return;
	visits[u] = true;
	renumbers[u] = renumber_cnt;
	real_numbers[renumber_cnt] = u;
	renumber_cnt++;
	for (int v : vertices[u])
		dfs(v);
}

void renumber() {
	std::priority_queue<ipair, std::vector<ipair>, std::greater<ipair>> Q;

	Q.emplace(0, rand() % V);
	while (!Q.empty()) {
		int u = Q.top().second;
		Q.pop();
		if (visits[u]) continue;
		visits[u] = true;
		prenumbers[u] = renumber_cnt++;
		for (int v : vertices[u]) {
			if (visits[v]) continue;
			degrees[v]++;
			Q.emplace(-degrees[v], v);
		}
	}
	for (int i = 0; i < V; i++)
		std::sort(vertices[i].begin(), vertices[i].end(),
			[](int x, int y) {
				return prenumbers[x] < prenumbers[y];
			}
		);
	// TODO: some vertices cannot be visited
	for (int i = 0; i < V; i++)
		visits[i] = false;
	renumber_cnt = 0;
	for (int i = 0; i < V; i++)
		dfs(i);
	for (int i = 0; i < E; i++) {
		renumbered_edges[i][0] = renumbers[edges[i][0]];
		renumbered_edges[i][1] = renumbers[edges[i][1]];
		renumbered_edges[i][2] = edges[i][2];
	}
}

void try_GA() {
	group = population();
	do {
		int num_crossover = NUM_CHILDREN / 4;
		for (int i = 0; i < num_crossover; i++) {
			int x = rand() % group.num_chrs;
			int y = rand() % group.num_chrs;
			group.children[i] = group.chrs[x]->crossover(group.chrs[y]);
		}
		for (int i = num_crossover; i < NUM_CHILDREN; i++) {
			int x = rand() % group.num_chrs;
			group.children[i] = group.chrs[x]->mutation();
		}
		group.replace();
	} while (get_time() - starts_at < 3);
	// } while (get_time() - starts_at + SPARE_TIME < V / 6.0);
}

void print_output() {
	chromosome *best = group.evals[0].chr;

	for (int i = 0; i < V; i++)
		answers[real_numbers[i]] = best->genes[i / 8] & (1 << (i % 8));
	for (int i = 0; i < V; i++)
		if (answers[i])
			printf("%d ", i);
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
