//  Windows
#ifdef _WIN32

#include <Windows.h>

double get_wall_time() {
	LARGE_INTEGER time, freq;
	if (!QueryPerformanceFrequency(&freq)) {
		//  Handle error
		return 0;
	}
	if (!QueryPerformanceCounter(&time)) {
		//  Handle error
		return 0;
	}
	return (double) time.QuadPart / freq.QuadPart;
}

double get_cpu_time() {
	FILETIME a, b, c, d;
	if (GetProcessTimes(GetCurrentProcess(), &a, &b, &c, &d) != 0) {
		//  Returns total user time.
		//  Can be tweaked to include kernel times as well.
		return (double) (d.dwLowDateTime | ((unsigned long long) d.dwHighDateTime << 32)) * 0.0000001;
	} else {
		//  Handle error
		return 0;
	}
}
//  Posix/Linux
#else
#include <time.h>
#include <sys/time.h>
double get_wall_time(){
	struct timeval time;
	if (gettimeofday(&time,NULL)){
		//  Handle error
		return 0;
	}
	return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
double get_cpu_time(){
	return (double)clock() / CLOCKS_PER_SEC;
}
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <stack>
#include <queue>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>

using namespace std;

class Knowledgebase {
private:
	//Define maximum amount of time for each query in seconds
	long maxQueryTime = 10;

	//Define meaning of argument
	typedef vector<string> argumentList;

	//Define meaning of predicate
	struct predicate {
		bool negated;
		string name;
		argumentList args;

		predicate negate() {
			bool s = !negated;
			string n = name;
			argumentList a = args;
			return predicate(s, n, a);
		}

		predicate(string &s) {
			const char *t = s.c_str();

			negated = t[0] == '~';

			unsigned int index = 0;
			while (t[index++] != '(') {};

			if (negated)
				name = string(t, 1, index - 2);
			else
				name = string(t, 0, index - 1);

			unsigned int length = (unsigned int) s.length();

			for (unsigned int i = index; i < length; i++) {
				if (t[i] == ',' || t[i] == ')') {
					string temp(t, index, i - index);
					index = i + 1;
					args.push_back(temp);
				}
			}
		}

		predicate(bool &s, string &n, argumentList &a) {
			negated = s;
			name = n;
			args = a;
		}

		bool operator==(const predicate &p) const {
			size_t signature1 = 17;
			signature1 = signature1 * 31 + hash<bool>()(this->negated);
			signature1 = signature1 * 31 + hash<string>()(this->name);
			for (int i = 0; i < this->args.size(); i++) {
				if (isVariable(this->args[i]))
					signature1 = signature1 * 31 + hash<string>()("x");
				else
					signature1 = signature1 * 31 + hash<string>()(this->args[i]);
			}

			size_t signature2 = 17;
			signature2 = signature2 * 31 + hash<bool>()(p.negated);
			signature2 = signature2 * 31 + hash<string>()(p.name);
			for (int i = 0; i < p.args.size(); i++) {
				if (isVariable(p.args[i]))
					signature2 = signature2 * 31 + hash<string>()("x");
				else
					signature2 = signature2 * 31 + hash<string>()(p.args[i]);
			}

			return signature1 == signature2;
		}
	};

	//Define meaning of sentence
	typedef vector<predicate> sentence;

	//Define hash function for predicates
	struct hash_predicate {
		size_t operator()(const predicate &p) const {
			size_t signature = 17;
			signature = signature * 31 + hash<bool>()(p.negated);
			signature = signature * 31 + hash<string>()(p.name);
			for (int j = 0; j < p.args.size(); j++) {
				signature = signature * 31 + hash<string>()(p.args[j]);
			}
			return signature;
		}
	};

	//Define hash function for sentences
	struct hash_sentence {
		size_t operator()(const sentence &s) const {
			//Variable rewriting is used to equalize the sentences
			//After rewriting the order of the variables in the predicates matters and not the variable names
			int variable_count = 0;
			unordered_map<string, int> visited_variables;

			size_t signature = 0;
			for (int i = 0; i < s.size(); i++) {
				size_t temp = 17;
				temp = temp * 31 + hash<bool>()(s[i].negated);
				temp = temp * 31 + hash<string>()(s[i].name);
				for (int j = 0; j < s[i].args.size(); j++) {
					if (isVariable(s[i].args[j])) {
						if (visited_variables.count(s[i].args[j]) == 0) {
							//If variable hasn't been visited before
							//Add to visited_variables
							variable_count++;
							visited_variables[s[i].args[j]] = variable_count;
						}
						temp = temp * 31 + hash<int>()(visited_variables[s[i].args[j]]);
					} else {
						temp = temp * 31 + hash<string>()(s[i].args[j]);
					}
				}
				signature ^= temp;
			}
			return signature;
		}
	};

	//Helper functions and classes
	static bool isVariable(const string &x) {
		return islower(x[0]);
	}

	static bool isLiteral(const predicate &p) {
		for (int i = 0; i < p.args.size(); i++) {
			if (isVariable(p.args[i]))
				return false;
		}

		return true;
	}

	class CNF {
	private:
		struct node {
			node *parent, *left, *right;
			string data;

			node() : parent(nullptr), left(nullptr), right(nullptr), data("") {}
		};

		static bool isOperator(string &s) {
			return (s == "(" || s == ")" || s == "~" ||
			        s == "&" || s == "|" || s == "=>");
		}

		static int operatorPrecedence(string &op) {
			/*if (op == "(" || op == ")")
				return 5;*/
			if (op == "~")
				return 4;
			if (op == "&")
				return 3;
			if (op == "|")
				return 2;
			if (op == "=>")
				return 1;
			return 0;
		}

		static node *deepCopy(node *root) {
			if (root == nullptr)
				return nullptr;

			node *temp = new node;
			temp->data = root->data;
			temp->left = deepCopy(root->left);
			if (temp->left != nullptr)
				temp->left->parent = temp;
			temp->right = deepCopy(root->right);
			if (temp->right != nullptr)
				temp->right->parent = temp;

			return temp;
		}

		string createExpressionString(node *root) {
			if (root == nullptr)
				return "";

			if (root->data == "~")
				return (root->data + " " + createExpressionString(root->left));
			else
				return (createExpressionString(root->left) + " " + root->data + " " +
				        createExpressionString(root->right));
		}

		static node *createExpressionTree(vector<string> expression) {
			stack<string> op;
			stack<node *> result;

			for (int i = 0; i < expression.size(); i++) {
				string temp = expression[i];
				if (isOperator(temp)) {
					if (op.empty() || (temp != ")" && op.top() == "(")) {
						op.push(temp);
					} else if (temp == "(") {
						op.push(temp);
					} else if (temp == ")") {
						while (op.top() != "(") {
							addToStack(op.top(), result);
							op.pop();
						}
						op.pop();
					} else if (operatorPrecedence(temp) > operatorPrecedence(op.top())) {
						op.push(temp);
					} else if (operatorPrecedence(temp) <= operatorPrecedence(op.top())) {
						while (!op.empty() && operatorPrecedence(temp) <= operatorPrecedence(op.top())) {
							addToStack(op.top(), result);
							op.pop();
						}
						op.push(temp);
					}
				} else {
					addToStack(temp, result);
				}
			}

			while (!op.empty()) {
				addToStack(op.top(), result);
				op.pop();
			}

			return result.top();
		}

		vector<string> createExpressionList(node *root) {
			vector<string> result, temp;
			if (root == nullptr)
				return result;

			temp = createExpressionList(root->left);
			if (!temp.empty())
				result.insert(result.end(), temp.begin(), temp.end());

			result.push_back(root->data);

			temp = createExpressionList(root->right);
			if (!temp.empty())
				result.insert(result.end(), temp.begin(), temp.end());

			return result;
		}

		static void deleteExpressionTree(node *root) {
			if (root == nullptr)
				return;

			deleteExpressionTree(root->left);
			deleteExpressionTree(root->right);
			delete root;
		}

		static vector<string> tokenize(string s) {
			s.erase(std::remove_if(s.begin(), s.end(), ::isspace), s.end());
			vector<string> result;

			const char *exp = s.c_str();
			unsigned length = (unsigned) s.length();

			for (unsigned i = 0; i < length; i++) {
				if (exp[i] == '(') {
					result.push_back("(");
				} else if (exp[i] == ')') {
					result.push_back(")");
				} else if (exp[i] == '~') {
					result.push_back("~");
				} else if (exp[i] == '&') {
					result.push_back("&");
				} else if (exp[i] == '|') {
					result.push_back("|");
				} else if (exp[i] == '=' && exp[i + 1] == '>') {
					result.push_back("=>");
					i++;
				} else {
					auto j = i + 1;
					while (exp[j] != ')') { j++; };
					string temp(exp, i, j - i + 1);
					result.push_back(temp);
					i = j;
				}
			}

			return result;
		}

		static void addToStack(string &op, stack<node *> &operandStack) {
			//create node for the operator or operand
			node *root = new node;
			root->data = op;

			if (isOperator(op)) {
				//if node is an operator process the operator based on unary or binary operation
				if (op == "~") {
					//get node for the operand
					node *operand = operandStack.top();
					operand->parent = root;
					//pop the operand stack
					operandStack.pop();
					//assign the operand to the left of the root operator
					root->left = operand;
					root->right = nullptr;
				} else {
					//get node for second operand
					node *operand2 = operandStack.top();
					operand2->parent = root;
					operandStack.pop();
					//get node for first operand
					node *operand1 = operandStack.top();
					operand1->parent = root;
					operandStack.pop();
					//assign the operand to the left and right of the root operator
					root->left = operand1;
					root->right = operand2;
				}
			}
			operandStack.push(root);
		}

		static node *negate(node *root) {
			if (root == nullptr)
				return nullptr;

			if (root->left == nullptr && root->right == nullptr) {
				//If leaf node is encountered i.e operand encountered
				if (root->data[0] != '~') {
					//premise is positive
					root->data = "~" + root->data;
				} else {
					//premise is negative
					root->data = root->data.substr(1);
				}
			} else {
				//operator encountered
				if (root->data == "~") {
					//If a "~" operator node is encountered
					node *t = root;
					if (root->parent == nullptr) {
						//if root doesnt have a parent then make child the root
						root = root->left;
						root->parent = nullptr;
					} else {
						if (root->parent->left == root) {
							//if "~" is on the left side of parent
							root->parent->left = root->left;
						} else {
							//if "~" is on the right side of parent
							root->parent->right = root->left;
						}
						//make the roots parent, the child's parent
						root->left->parent = root->parent;
						root = root->left;
						//delete the "~" node
					}
					delete t;
					return root;
				} else if (root->data == "&") {
					root->data = "|";
				} else if (root->data == "|") {
					root->data = "&";
				}
			}

			root->left = negate(root->left);
			root->right = negate(root->right);

			return root;
		}

		static void removeImplications(node *root) {
			if (root == nullptr)
				return;
			removeImplications(root->left);
			removeImplications(root->right);
			if (root->data == "=>") {
				//A => B ------> ~A | B
				negate(root->left);
				root->data = "|";
			}
		}

		static node *resolveNegations(node *root) {
			if (root == nullptr)
				return nullptr;

			while (root->data == "~") {
				root = negate(root); //removes the "~" sign only
				root = negate(root); //performs the actual negation
			}

			root->left = resolveNegations(root->left);
			root->right = resolveNegations(root->right);

			return root;
		}

		static node *distribute(node *parent, node *child) {
			node *grandparent = parent->parent;
			node *leftBranch1, *leftBranch2, *rightBranch1, *rightBranch2;
			bool isParentLeftOfGrandparent = (grandparent != nullptr) ? (grandparent->left == parent) : false;

			if (parent->left == child) {
				//if child is the left child of the parent
				leftBranch1 = child->left;
				leftBranch2 = child->right;
				rightBranch1 = parent->right;
				rightBranch2 = deepCopy(parent->right);
			} else {
				//If child is the right child of the parent
				leftBranch1 = parent->left;
				leftBranch2 = deepCopy(parent->left);
				rightBranch1 = child->left;
				rightBranch2 = child->right;
			}

			//delete parent
			delete parent;

			//create and initialize new left and right parent
			node *leftNode = new node;
			leftNode->data = "|";
			leftNode->left = leftBranch1;
			leftNode->right = rightBranch1;
			leftNode->parent = child;

			node *rightNode = new node;
			rightNode->data = "|";
			rightNode->left = leftBranch2;
			rightNode->right = rightBranch2;
			rightNode->parent = child;

			//Make left node parent of *branch1
			leftBranch1->parent = leftNode;
			rightBranch1->parent = leftNode;

			//Make right node parent of *branch2
			leftBranch2->parent = rightNode;
			rightBranch2->parent = rightNode;

			//Make child point to new left and right nodes
			child->left = leftNode;
			child->right = rightNode;

			//Make child grandparents child
			child->parent = grandparent;
			if (grandparent != nullptr) {
				//If child has a grandparent
				if (isParentLeftOfGrandparent) {
					grandparent->left = child;
				} else {
					grandparent->right = child;
				}
			}
			return child;
		}

		static node *distributeOrOverAnd(node *root) {
			if (root == nullptr)
				return nullptr;

			//If root is an operand, return without modification
			if (!isOperator(root->data))
				return root;

			//If root is "|"...
			if (root->data == "|") {
				bool distributed = false;
				//...and its left child is "&"
				if (root->left->data == "&") {
					//distribute | over & on the left child
					root = distribute(root, root->left);
					distributed = true;
				}

				//...and its right child is "&"
				if (root->right->data == "&") {
					//distribute | over & on the right child
					root = distribute(root, root->right);
					distributed = true;
				}

				//If distribution has taken place
				if (distributed) {
					if (root->parent == nullptr) {
						//If root is the absolute root of the tree, test children for possible distributions
						return distributeOrOverAnd(root);
					} else {
						//Bubble up the control flow
						return root;
					}
				}
			}

			//Store original left and right children
			node *left = root->left;
			node *right = root->right;
			/*
			 *  If after calling distributeOrOverAnd on the left or right subtree,
			 *  the left or right subtree changes then we need to run again for possible distributions in the children
			 */
			root->left = distributeOrOverAnd(root->left);
			if (left != root->left) {
				//If left child changes then return the result of rerunning over root
				return distributeOrOverAnd(root);
			}
			root->right = distributeOrOverAnd(root->right);
			if (right != root->right) {
				//If right child changes then return the result of rerunning over root
				return distributeOrOverAnd(root);
			}

			return root;
		}

		static vector<node *> splitSentenceOverAnd(node *root) {
			vector<node *> result, temp;
			if (root == nullptr)
				return result;

			if (root->data != "&") {
				result.push_back(root);
				return result;
			}

			temp = splitSentenceOverAnd(root->left);
			if (!temp.empty()) {
				result.insert(result.end(), temp.begin(), temp.end());
			}
			temp = splitSentenceOverAnd(root->right);
			if (!temp.empty()) {
				result.insert(result.end(), temp.begin(), temp.end());
			}

			return result;
		}

		static sentence createCNFSentence(node *root) {
			sentence result, temp;
			if (root->left == nullptr && root->right == nullptr) {
				predicate t = predicate(root->data);
				result.push_back(t);
				return result;
			}

			temp = createCNFSentence(root->left);
			if (!temp.empty())
				result.insert(result.end(), temp.begin(), temp.end());

			temp = createCNFSentence(root->right);
			if (!temp.empty())
				result.insert(result.end(), temp.begin(), temp.end());

			return result;
		}

	public:
		static sentence negateCNFSentence(sentence s) {
			for (int i = 0; i < s.size(); i++) {
				s[i] = s[i].negate();
			}
			return s;
		}

		static sentence &factorize(sentence &s) {
			//Removes duplicate predicates
			//Set of visited predicates
			unordered_set<predicate, hash_predicate> visited;
			//Vector of predicate indices to remove
			vector<int> predicate_to_remove;
			for (int i = 0; i < s.size(); i++) {
				if (visited.count(s[i]) == 0) {
					//If predicate hasn't been encountered before
					//Mark as visited
					visited.insert(s[i]);
				} else {
					//Else predicate has been encountered before
					//Mark for removal
					predicate_to_remove.push_back(i);
				}
			}
			//Remove marked indices from sentence
			for (int i = 0; i < predicate_to_remove.size(); i++) {
				s.erase(s.begin() + predicate_to_remove[i]);
			}

			return s;
		}

		static vector<sentence> convertToCNFSentences(string &s) {
			//Create expression tree from tokenized string
			node *expressionRoot = createExpressionTree(tokenize(s));
			//cout << "Tokenized Expression \t\t:\t " << createExpressionString(expressionRoot) << endl;

			//Remove implications
			removeImplications(expressionRoot);
			//cout << "Expression without Implications :\t " << createExpressionString(expressionRoot) << endl;

			//Resolve negations
			expressionRoot = resolveNegations(expressionRoot);
			//cout << "Expression without Negations \t:\t " << createExpressionString(expressionRoot) << endl;

			//Distribute | over &
			expressionRoot = distributeOrOverAnd(expressionRoot);
			//cout << "After Distribution of | over & \t:\t " << createExpressionString(expressionRoot) << endl;

			//Split sentences over &
			vector<node *> sentences = splitSentenceOverAnd(expressionRoot);
			//cout << "New Sentences in CNF" << endl;

			sentence temp;
			vector<sentence> result;
			for (int i = 0; i < sentences.size(); i++) {
				//Convert to CNF sentence
				temp = createCNFSentence(sentences[i]);
				//Factorize the CNF sentence
				//Removes redundant predicates
				//A(x) | A(x) becomes A(x)
				temp = factorize(temp);
				//Add to result
				result.push_back(temp);
			}

			//Delete the expression tree
			deleteExpressionTree(expressionRoot);

			return result;
		}
	};

	//Define Database class
	class Database {
	private:
		struct row {
			vector<pair<unsigned int, unsigned int>> positive_literals;
			vector<pair<unsigned int, unsigned int>> negative_literals;
			vector<pair<unsigned int, unsigned int>> positive_sentences;
			vector<pair<unsigned int, unsigned int>> negative_sentences;
		};

		unordered_map<char, unsigned int> variables; //Used for standardization

		vector<sentence> data;  //Contains all sentences
		unordered_map<string, row> index; //Contains index about all sentences

		sentence standardizeSentence(sentence &s) {
			unordered_set<char> current_variables;
			//Iterate through all predicates of the sentence
			for (int i = 0; i < s.size(); i++) {
				argumentList &args = s[i].args;
				//Iterate through all arguments of the predicate
				for (int j = 0; j < args.size(); j++) {
					if (isVariable(args[j])) {
						//If argument is a variable then standardize
						char var = args[j][0];
						if (current_variables.count(var) == 0) {
							//If variable hasn't been seen in the sentence previously
							if (variables.count(var) == 0) {
								//If variable hasn't been seen in the KB previously
								//Initialize the variable's unique count
								variables[var] = 1;

							} else {
								//If variable has been seen in the KB previously
								//Increment the variable's unique count
								variables[var]++;
							}
							//insert the variable in the current variable set
							current_variables.insert(var);
						}
						//Append unique id
						args[j] = var + to_string(variables[var]);
					}
				}
			}

			return s;
		}

	public:
		Database copy() {
			Database t;
			t.variables = this->variables;
			t.data = this->data;
			t.index = this->index;
			return t;
		}

		void store(sentence &s) {
			//Standardize the CNF sentence
			s = standardizeSentence(s);
			data.push_back(s);
			//Generate location index of the sentence in the data vector
			unsigned long long int loc = data.size() - 1;
			//Loop through all predicates in the sentence
			for (int i = 0; i < s.size(); i++) {
				//Index the predicate
				if (isLiteral(s[i])) {
					if (!s[i].negated) {
						index[s[i].name].positive_literals.push_back(pair<unsigned int, unsigned int>(loc, i));
					} else {
						index[s[i].name].negative_literals.push_back(pair<unsigned int, unsigned int>(loc, i));
					}
				} else {
					if (!s[i].negated) {
						index[s[i].name].positive_sentences.push_back(pair<unsigned int, unsigned int>(loc, i));
					} else {
						index[s[i].name].negative_sentences.push_back(pair<unsigned int, unsigned int>(loc, i));
					}
				}
			}
		}

		vector<pair<sentence, unsigned int>> fetch(predicate &p) {
			vector<pair<sentence, unsigned int>> result;
			vector<pair<unsigned int, unsigned int>> literalIndex;
			vector<pair<unsigned int, unsigned int>> sentenceIndex;

			if (!p.negated) {
				literalIndex = index[p.name].positive_literals;
				sentenceIndex = index[p.name].positive_sentences;
			} else {
				literalIndex = index[p.name].negative_literals;
				sentenceIndex = index[p.name].negative_sentences;
			}
			//Fetch all literals and add to result
			for (int i = 0; i < literalIndex.size(); i++) {
				result.push_back(pair<sentence, unsigned int>(data[literalIndex[i].first], literalIndex[i].second));
			}
			//Fetch all non literals and add to result
			for (int i = 0; i < sentenceIndex.size(); i++) {
				result.push_back(pair<sentence, unsigned int>(data[sentenceIndex[i].first], sentenceIndex[i].second));
			}

			return result;
		}

	};

	//Create a global database for sentences in the knowledge base
	Database DB;

	//Substitutes argumentList with substitution list theta
	argumentList &substitute(argumentList &x, unordered_map<string, string> &theta) {
		//Loop through all arguments
		for (int i = 0; i < x.size(); i++) {
			//Check if argument is substitutable using theta
			while (theta.count(x[i]) > 0)
				x[i] = theta[x[i]];
		}
		return x;
	}

	//Unifies argumentList x with argumentList y and creates a substitution list theta
	bool unify(argumentList &x, argumentList &y, unordered_map<string, string> &theta) {
		if (x.size() != y.size())
			return false;
		for (int i = 0; i < x.size(); i++) {
			if (x[i] != y[i]) {
				if (isVariable(x[i])) {
					//If x[i] is a variable
					theta[x[i]] = y[i];
					x = substitute(x,theta);
					y = substitute(y,theta);
				} else if (isVariable(y[i])) {
					//If y[i] is a variable
					theta[y[i]] = x[i];
					x = substitute(x,theta);
					y = substitute(y,theta);
				} else {
					//If x[i] and y[i] both are constants
					return false;
				}
			}
		}
		return true;
	}

public:
	void tell(string &fact) {
		//Loop through all CNF sentences and insert into KB
		vector<sentence> sentences = CNF::convertToCNFSentences(fact);
		for (int i = 0; i < sentences.size(); i++) {
			DB.store(sentences[i]);
		}
	}

	bool ask(string &query) {
		//Calculate finish time quota
		double finishTime = get_wall_time() + maxQueryTime;
		//Create query into CNF Sentence
		sentence alpha = CNF::convertToCNFSentences(query)[0];
		//Negate alpha
		sentence notAlpha = CNF::negateCNFSentence(alpha);
		//Clone the knowledge base's data
		Database KB = DB.copy();
		//Store alpha in KB (KB ^ ~alpha = unsatisfiable)
		KB.store(notAlpha);

		queue<sentence> Frontier;
		unordered_set<sentence, hash_sentence> LoopDetector;    //Prevents duplicate sentences in KB

		Frontier.push(notAlpha);

		while (!Frontier.empty()) {
			//Choose shallowest node from frontier
			sentence currentSentence = Frontier.front();
			//Remove the node from frontier
			Frontier.pop();

			for (int i = 0; i < currentSentence.size(); i++) {

				//create a resolver predicate by negating the current predicate
				predicate resolver = currentSentence[i].negate();
				//Get resolvableSentences for each predicate in the currentSentence
				vector<pair<sentence, unsigned int>> resolvableSentences = KB.fetch(resolver);
				for (int j = 0; j < resolvableSentences.size(); j++) {
					//Resolve each sentence in the resolvableSentences
					//Create substitution list
					unordered_map<string, string> theta;
					//Find unifiable predicate in the sentence
					if (resolvableSentences[j].first[resolvableSentences[j].second].name == currentSentence[i].name &&
					    resolvableSentences[j].first[resolvableSentences[j].second].negated !=
					    currentSentence[i].negated) {

						argumentList x = currentSentence[i].args;
						argumentList y = resolvableSentences[j].first[resolvableSentences[j].second].args;

						if (unify(x, y, theta)) {
							//If unifiable, use substitution list theta to unify all the predicates in both sentences
							sentence t1 = currentSentence;
							sentence t2 = resolvableSentences[j].first;
							for (int k = 0; k < t1.size(); k++)
								t1[k].args = substitute(t1[k].args, theta);
							for (int k = 0; k < t2.size(); k++)
								t2[k].args = substitute(t2[k].args, theta);

							//Resolve sentences t1 and t2;

							//Remove the resolving predicates
							t1.erase(t1.begin() + i);
							t2.erase(t2.begin() + resolvableSentences[j].second);

							//Merge the two resolved sentences to obtain the resolvent
							sentence resolvent;
							resolvent.insert(resolvent.end(), t1.begin(), t1.end());
							resolvent.insert(resolvent.end(), t2.begin(), t2.end());

							//Factorize the resolvent to remove duplicate predicates
							resolvent = CNF::factorize(resolvent);

							//If resolvent is empty then alpha is true
							if (resolvent.empty()) {
								DB.store(alpha);
								return true;
							}

							if (LoopDetector.count(resolvent) == 0) {
								//If resolvent hasn't been encountered before then add to Frontier and store in KB
								KB.store(resolvent);
								Frontier.push(resolvent);
								LoopDetector.insert(resolvent);
							}
						}
					}
				}
				//Check if allotted time ran out
				if (get_wall_time() > finishTime) {
					return false;
				}
			}
		}
		return false;
	}

};

int main() {
	int NQ = 0, NS = 0;
	Knowledgebase KB;
	string tempString;
	vector<string> query;

	ifstream InputFile("input.txt");
	if (InputFile.is_open()) {
		getline(InputFile, tempString);
		NQ = stoi(tempString);

		query = vector<string>((unsigned) NQ);
		for (int i = 0; i < NQ; i++) {
			getline(InputFile, tempString);
			query[i] = tempString;
		}

		getline(InputFile, tempString);
		NS = stoi(tempString);

		for (int i = 0; i < NS; i++) {
			getline(InputFile, tempString);
			KB.tell(tempString);
		}
		InputFile.close();
	} else {
		cout << "Input file failed to load" << endl;
	}

	string output = "";
	for (int i = 0; i < NQ; i++) {
		bool t = KB.ask(query[i]);
		output += t ? "TRUE" : "FALSE";
		output += "\n";
	}

	ofstream OutputFile("output.txt");
	if (OutputFile.is_open()) {
		OutputFile << output;
		cout << "Output Success" << endl;
		OutputFile.close();
	} else {
		cout << "Output file failed to load" << endl;
	}

	return 0;
}
