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
std::vector<int> vertices[MAX_V];  // only used while renumbering

int hash_const;  // (2^V - 1) % MOD
int renumber_cnt;
bool visits[MAX_V];
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
	unsigned char genes[MAX_L];

	chromosome(bool is_random = false) {
		if (is_random)
			for (int i = 0; i < L; i++)
				genes[i] = rand() % 256;
	}

	chromosome(chromosome *other) {
		for (int i = 0; i < L; i++)
			genes[i] = other->genes[i];
	}

	// copy interval [left, (right + 7) / 8 * 8) from other to this
	void get_interval(chromosome *other, int left, int right, bool flip) {
		int pos1 = left / 8, pos2 = left % 8;
		genes[pos1] = (genes[pos1] & ((1 << pos2) - 1))
			    | ((flip ? ~other->genes[pos1] : other->genes[pos1]) & (255 << pos2));
		for (int i = pos1 + 1; i <= (right - 1) / 8; i++)
			genes[i] = (flip ? ~other->genes[i] : other->genes[i]);
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
		bool flip = rand() % 2;
		for (int i = 0; i <= CUTTING_POINT; i++)
			child->get_interval((i % 2) ? other : this, cp[i], cp[i + 1], (i % 2) && flip);
		return child;
	}

	chromosome *mutation() {
		chromosome *child = new chromosome(this);
		int idx = rand() % V;
		child->genes[idx / 8] ^= 1 << (idx % 8);
		return child;
	}

	int hash() const {
		int result = 0;
		for (int i = 0; i < L - 1; i++)
			result = (256LL * result + genes[i]) % MOD;
		result = (256LL * result + (genes[L - 1] & ((1 << ((V - 1) % 8 + 1)) - 1))) % MOD;
		return std::min((int)result, get_complement((int)result));
	}

	int evaluate() const {
		for (int i = 0; i < V; i++)
			visits[i] = genes[i / 8] >> (i % 8) & 1;
		int score = 0;
		for (int i = 0; i < E; i++)
			if (visits[edges[i][0]] != visits[edges[i][1]])
			    score += edges[i][2];
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
} group;

void get_input() {
	// get input
	int u, v, w;
	if (scanf("%d %d", &V, &E) != 2) exit(errno);
	L = (V + 7) / 8;
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
		renumbers[u] = renumber_cnt++;
		for (int v : vertices[u]) {
			if (visits[v]) continue;
			degrees[v]++;
			Q.emplace(-degrees[v], v);
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
	for (int i = 0; i < V; i++)
		visits[i] = false;
	renumber_cnt = 0;
	for (int i = 0; i < V; i++)
		dfs(i);
	for (int i = 0; i < E; i++) {
		edges[i][0] = renumbers[edges[i][0]];
		edges[i][1] = renumbers[edges[i][1]];
	}
}

void try_GA() {
	group = population();
	// int cnt = 0;
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
		// cnt++;
		// if (cnt % 100 == 0)
		// 	fprintf(stderr, "%d %d %lf\n", cnt, group.evals[0].score, get_time() - starts_at);
	} while (get_time() - starts_at + SPARE_TIME < V / 6.0);
}

void print_output() {
	chromosome *best = group.evals[0].chr;

	for (int i = 0; i < V; i++)
		visits[real_numbers[i]] = (best->genes[i / 8] >> (i % 8)) & 1;
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
