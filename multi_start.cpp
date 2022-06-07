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

#define MAX_POPULATION 8175616

int V, E;
int edges[MAX_E][3];
std::vector<int> vertices[MAX_V];
std::vector<std::pair<int, int>> infos[MAX_V];
std::priority_queue<std::pair<int, int>> Q;

int renumber_cnt;
uint8_t visits[MAX_V];
int degrees[MAX_V], renumbers[MAX_V], real_numbers[MAX_V];

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
	chromosome *chr;

	evaluation() = default;
	evaluation(chromosome *chr_) : chr(chr_) {
		score = chr->evaluate();
	}

	bool operator<(const evaluation &other) {
		return score > other.score;
	}
};

class population {
public:
	chromosome *chrs[MAX_POPULATION];
	evaluation evals[MAX_POPULATION];

	population() {
		for (int i = 0; i < MAX_POPULATION; i++) {
			chrs[i] = new chromosome(true);
			chrs[i]->local_opt();
			evals[i] = evaluation(chrs[i]);
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

int main() {
	// srand, rand is fast, we do not need true-randomness
	srand(time(NULL));

	get_input();
	renumber();
	group = new population();
	int score = 0;
	for (int i = 0; i < MAX_POPULATION; i++)
		score = std::max(score, group->evals[i].score);
	fprintf(stderr, "%d\n", score);
	return 0;
}